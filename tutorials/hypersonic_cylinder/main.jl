# main.jl - Entry point for hypersonic cylinder simulation
#
# Description:
#   Main driver script for running a hypersonic flow stability simulation
#   over a cylinder geometry. This script orchestrates the full simulation
#   pipeline:
#     1. Activates the project environment and loads the solver module
#     2. Loads user-defined input parameters from input_file.jl
#     3. Initializes the chemistry model (if not restarting)
#     4. Sets up the computational mesh and initial conditions
#     5. Runs the time-marching simulation (multiple stages)
#     6. Performs linearization, eigenvalue, and transient growth analysis
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
solution_save = nothing
solution = INITIALIZATION(solution, solution_save, chemistry)
PLOT_INITIAL_SET_UP(solution)

## Run simulation (coarse grid lower speed for convergence)
solution["restart"]        = false      # Start from scratch
solution["time_integration"]["N_iter"]   = 1500
solution["freestream"]["u"]   = 2000       # [m/s]
solution["freestream"]["Re"]  = 2000       # Reynolds number (based on tip radius)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

## Run simulation (coarse grid increase speed)
# Need to increase gradually freestream speed, because if not strong shock can be
# generated after the main shock and can break numerics
solution["restart"]        = true      # Restart from previous case
solution["shock"]["interpolate"]              = "1st" # Set first order interpolation for increased robustness in strong shock transient
solution["time_integration"]["N_iter"]   = 3500
solution["freestream"]["u"]   = 5000     # [m/s]
solution["freestream"]["Re"]  = 2000     # Reynolds number (based on tip radius)
solution = INITIALIZATION(solution, solution_save, chemistry) # Restart from previous case

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

## Run simulation (remesh to finer grid & higher Reynolds)
solution["restart"]        = true      # Restart from previous case
solution["shock"]["interpolate"]              = "2nd"
solution["time_integration"]["N_iter"]   = 1000 # Increase more iterations in case want better converged case
solution["freestream"]["u"]   = 5000      # [m/s]
solution["freestream"]["Re"]  = 20000     # Reynolds number (based on tip radius)
solution["mesh"]["Nchi"] = 400  # Points along chi (wall boundary)
solution["mesh"]["Neta"] = 80   # Points along eta (wall-normal direction)
if !haskey(solution["curvilinear_mapping"], "refinement_stagnation")
    solution["curvilinear_mapping"]["refinement_stagnation"] = Dict{String, Any}()
end
solution["curvilinear_mapping"]["refinement_stagnation"]["state"]        = true
solution["curvilinear_mapping"]["refinement_stagnation"]["BL_thickness"] = 0.2 # Thickness region that is refined
solution["curvilinear_mapping"]["refinement_stagnation"]["intensity"]    = 0.95 # 0 - not refined, 1 - infinitely refined
if !haskey(solution["curvilinear_mapping"], "refinement_wall")
    solution["curvilinear_mapping"]["refinement_wall"] = Dict{String, Any}()
end
solution["curvilinear_mapping"]["refinement_wall"]["state"]              = true
solution["curvilinear_mapping"]["refinement_wall"]["BL_thickness"]       = 0.1 # Thickness region that is refined
solution["curvilinear_mapping"]["refinement_wall"]["intensity"]          = 0.995 # 0 - not refined, 1 - infinitely refined
solution = INITIALIZATION(solution, solution_save, chemistry) # Restart from previous case

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

## Visualize density
plot_x_coords = solution["mesh"]["x"]
plot_y_coords = solution["mesh"]["y"]
plot_c_data = solution["var"]["rho"][2:end-1, 2:end-1] ./ solution["freestream"]["rho_0"]
(fig1, ax1) = curvilinear_heatmap(plot_x_coords, plot_y_coords, plot_c_data,
                    colorbar_label=L"$\rho/\rho_\infty$",
                    xlabel_str="x/R", ylabel_str="y/R",
                    fig_size=auto_figure_size(plot_x_coords, plot_y_coords))
fig1.tight_layout()
save_or_display_figure(fig1, name="cylinder_density")

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
save_or_display_figure(fig2, name="cylinder_vorticity")

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
PLOT_MODES(freestream_disturbances, L, solution, chemistry, V[:, mode], T_plot) # Plot modes

## Compute transient growth downstream flow
T_TGD = [0.5; 1.0; 1.5] # Time to find optimal transient growth
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
time_optimization_index = 2
mode = 1 # Mode to visualize
T_plot = 1.50 # Time in which disturbance is plotted
freestream_disturbances = false
PLOT_MODES(freestream_disturbances, L, solution, chemistry, V_TGD[:, mode, time_optimization_index], T_plot) # Plot modes

## Linearize A_ (extended system with freestream perturbations)
N_l = 40
w_infty = 2*pi*collect(range(0, N_l, length=N_l+1)) # Streamwise frequencies allowed of freestream disturbances

(solution, L_) = LINEARIZE_L_(L, solution, chemistry, w_infty)

## Compute freestream receptivity optimal modes
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
mode = 2 # Mode to visualize
T_plot = 5 # Time in which disturbance is plotted
freestream_disturbances = true
PLOT_MODES(freestream_disturbances, L_, solution, chemistry, V_TGF[:, mode, time_optimization_index], T_plot, w_infty) # Plot modes
