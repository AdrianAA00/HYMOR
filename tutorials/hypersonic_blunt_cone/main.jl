# main.jl - Entry point for 3D axisymmetric blunted cone simulation
#
# Description:
#   Main driver script for running a hypersonic flow stability simulation
#   over a 3D axisymmetric blunted cone geometry. This script orchestrates
#   the full simulation pipeline:
#     1. Activates the project environment and loads the solver module
#     2. Loads user-defined input parameters from input_file.jl
#     3. Initializes the chemistry model (if not restarting)
#     4. Sets up the computational mesh and initial conditions
#     5. Runs the time-marching simulation (coarse + fine grid)
#     6. Performs linearization, eigenvalue, transient growth, and
#        freestream receptivity analysis
#     7. Generates post-processing plots
#
# Inputs:
#   None (all parameters are loaded from input_file.jl)
#
# Outputs:
#   solution     - Dict : Final solution state containing flow fields,
#                          mesh data, and solver metadata
#   solution_old - Dict : Solution from the previous time step
#   chemistry    - Dict : Chemistry model data (species, reactions, etc.)
#
# Usage:
#   Run this script from the test case directory:
#       julia main.jl
#
# Part of: Hypersonics Stability Solver - Test Cases Module

## Activate project environment and load solver module
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using HYMOR
using Printf
using LaTeXStrings
using Plots
using LinearAlgebra

## Solver root directory (for data file paths)
solver_dir = joinpath(@__DIR__, "..", "..") * "/"
solution = Dict{String, Any}()
solution["solver_dir"] = solver_dir

## Configure plot output directory (all figures saved here)
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output"))

## Load input configuration
filename = joinpath(@__DIR__, "input_file.jl")
LOAD_INPUT_VARIABLES(filename)

## Load chemistry model
if !solution["restart"]
    chemistry = SET_CHEMISTRY(solution)
end

## Plot initial setup
solution = INITIALIZATION(solution, solution_save, chemistry)
PLOT_INITIAL_SET_UP(solution)

## Run simulation (coarse grid for fast convergence)
solution["restart"]        = false     # Start from scratch
solution["time_integration"]["N_iter"]   = 2000
solution["freestream"]["Re"]  = 10000     # Reynolds number
solution["mesh"]["Nchi"] = 200  # Points along chi (wall boundary)
solution["mesh"]["Neta"] = 40   # Points along eta (wall-normal direction)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

## Run simulation (fine grid)
solution["restart"]        = true     # Start from scratch
solution["time_integration"]["N_iter"]   = 1000
solution["freestream"]["Re"]  = 100000     # Reynolds number
solution["mesh"]["Nchi"] = 800  # Points along chi (wall boundary)
solution["mesh"]["Neta"] = 100   # Points along eta (wall-normal direction)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

# ## Run simulation (finest grid), if run with this takes a while
# solution["restart"]        = true     # Start from scratch
# solution["time_integration"]["N_iter"]   = 1000
# solution["freestream"]["Re"]  = 100000     # Reynolds number
# solution["mesh"]["Nchi"] = 1600  # Points along chi (wall boundary)
# solution["mesh"]["Neta"] = 200   # Points along eta (wall-normal direction)
# solution = INITIALIZATION(solution, solution_save, chemistry)

# solution_base = deepcopy(solution)
# t_start = time()
# solution = RUN_SIMULATION(solution, solution_base, chemistry)
# t_elapsed = time() - t_start
# println("Elapsed time: $t_elapsed seconds")
# solution_save = deepcopy(solution)

## Visualize density
plot_x_coords = solution["mesh"]["x"]
plot_y_coords = solution["mesh"]["y"]
plot_c_data = solution["var"]["rho"][2:end-1, 2:end-1] ./ solution["freestream"]["rho_0"]
(fig1, ax1) = curvilinear_heatmap(plot_x_coords, plot_y_coords, plot_c_data,
                    colorbar_label=L"$\rho/\rho_\infty$",
                    xlabel_str="x/R", ylabel_str="y/R",
                    fig_size=auto_figure_size(plot_x_coords, plot_y_coords))
fig1.tight_layout()
save_or_display_figure(fig1, name="capsule_density")

## Visualize vorticity
plot_x_coords = solution["mesh"]["x"][2:end-1, 2:end-1]
plot_y_coords = solution["mesh"]["y"][2:end-1, 2:end-1]
u_vel_int = solution["var"]["rho_u"][2:end-1, 2:end-1] ./ solution["var"]["rho"][2:end-1, 2:end-1]
v_vel_int = solution["var"]["rho_v"][2:end-1, 2:end-1] ./ solution["var"]["rho"][2:end-1, 2:end-1]
(_, d_u_dy) = DERIVATIVE(u_vel_int, solution)
(d_v_dx, _) = DERIVATIVE(v_vel_int, solution)
plot_c_data = d_v_dx - d_u_dy
plot_c_data = plot_c_data / solution["freestream"]["U"] * solution["curvilinear_mapping"]["R"]
(fig2, ax2) = curvilinear_heatmap(plot_x_coords, plot_y_coords, plot_c_data,
                    colorbar_label=L"$\omega R/U_\infty$",
                    xlabel_str="x/R", ylabel_str="y/R",
                    fig_size=auto_figure_size(plot_x_coords, plot_y_coords))
fig2.tight_layout()
save_or_display_figure(fig2, name="capsule_vorticity")

## Linearize A
(solution, L) = LINEARIZE_L(solution, chemistry)

## Compute eigenvalues of stability analysis
n_modes = 10 # search 10 most unstable eigenmodes
(V, D) = EIGENVALUES(L, solution, n_modes) # V stores eigenmodes, and D the eigenvalues

## Visualize eigenmodes global modal analysis
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output", "modal"))
mode = 3 # Mode to visualize
T_plot = 0 # Time in which disturbance is plotted
freestream_disturbances = false
solution["running_plot"]["scaling_range"] = 1/2 # Saturate to see field
PLOT_MODES(freestream_disturbances, L, solution, chemistry, V[:, mode], T_plot) # Plot modes

## Compute eigenvalues of transient growth downstream flow
T_TGD = [3] # Time to find optimal transient growth
n_modes = 5 # Number of modes to find in each iteration
V_TGD = zeros(4 * solution["mesh"]["Nchi"] * solution["mesh"]["Neta"], n_modes, size(T_TGD, 1)) # Store modes
D_TGD = zeros(n_modes, n_modes, length(T_TGD)) # Store transient growth

T_opt_TGD = zeros(length(T_TGD), 1)
for i in 1:length(T_TGD)
    (V_i, D_i, T_opt_TGD[i, 1]) = TRANSIENT_GROWTH_DOWNSTREAM(L, solution, n_modes, T_TGD[i, 1])
    V_TGD[:, :, i] = V_i
    D_TGD[:, :, i] = diagm(D_i)
end

## Visualize eigenmodes  global non-modal analysis
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output", "non_modal"))
time_optimization_index = 1
mode = 2 # Mode to visualize
T_plot = 0.0 # Time in which disturbance is plotted
freestream_disturbances = false
solution["running_plot"]["scaling_range"] = 1/2 # Saturate to see field
PLOT_MODES(freestream_disturbances, L, solution, chemistry, V_TGD[:, mode, time_optimization_index], T_plot) # Plot modes

## Visualize selected non-modal
T_f = 5 # Final time of integration to see energy growth
mode = 2 # Mode to evolve in time
freestream_disturbances = false
get_amplification_only = true
max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances, V_TGD[:, mode, time_optimization_index], D_TGD[mode, mode, time_optimization_index], T_opt_TGD[time_optimization_index, 1], L, solution, chemistry, T_f, get_amplification_only)

println("max_ref(E) = " * string(max_gain["non_temporal"]["all"]))
println("max_ref(E_p) = " * string(max_gain["non_temporal"]["acoustic"]))
println("max_ref(E_S) = " * string(max_gain["non_temporal"]["entropic"]))
println("max_ref(E_k) = " * string(max_gain["non_temporal"]["kinetic"]))

## Linearize A_ (extended system with freestream perturbations)
N_l = 40
w_infty = 2*pi*collect(range(0, N_l, length=N_l+1)) # Streamwise frequencies allowed of freestream disturbances
(solution, L_) = LINEARIZE_L_(L, solution, chemistry, w_infty)

## Compute optimal freestream receptivity modes
T_TGF = [5] # Time to find optimal transient growth
n_modes = 5
V_TGF = zeros(ComplexF64, 4 * solution["mesh"]["Nchi"] * solution["mesh"]["Neta"] + 4 * solution["mesh"]["Nchi"] * (N_l + 1), n_modes, length(T_TGF))
D_TGF = zeros(ComplexF64, n_modes, n_modes, length(T_TGF))

T_opt_TGF = zeros(length(T_TGF), 1)
for i in 1:length(T_TGF)
    (V_i, D_i, T_opt_TGF[i, 1]) = FREESTREAM_RECEPTIVITY(L_, solution, n_modes, T_TGF[i, 1], w_infty)
    V_TGF[:, :, i] = V_i
    D_TGF[:, :, i] = diagm(D_i)
end

## Visualize freestream receptivity modes
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output", "freestream_receptivity"))
time_optimization_index = 1
mode = 3 # Mode to visualize
T_plot = 10 # Time in which disturbance is plotted
freestream_disturbances = true
solution["running_plot"]["scaling_range"] = 1/10 # Saturate to see field
PLOT_MODES(freestream_disturbances, L_, solution, chemistry, V_TGF[:, mode, time_optimization_index], T_plot, w_infty) # Plot modes

## Visualize selected freestream gains
T_f = 10 # Final time of integration to see energy growth
mode = 3 # Mode to evolve in time
freestream_disturbances = true
get_amplification_only = true
max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances, V_TGF[:, mode, time_optimization_index], D_TGF[mode, mode, time_optimization_index], T_opt_TGF[time_optimization_index, 1], L_, solution, chemistry, T_f, get_amplification_only, w_infty)

println("max_ref(E) = " * string(max_gain["non_temporal"]["all"]))
println("max_ref(E_p) = " * string(max_gain["non_temporal"]["acoustic"]))
println("max_ref(E_S) = " * string(max_gain["non_temporal"]["entropic"]))
println("max_ref(E_k) = " * string(max_gain["non_temporal"]["kinetic"]))
