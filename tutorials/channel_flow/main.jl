# main.jl - Entry point for channel flow stability simulation
#
# Description:
#   Main driver script for running a stability analysis on a channel flow.
#   This script orchestrates the full simulation pipeline:
#     1. Activates the project environment and loads the solver module
#     2. Loads user-defined input parameters from input_file.jl
#     3. Initializes the chemistry model
#     4. Sets up the computational mesh and initial conditions
#     5. Performs linearization and eigenvalue analysis
#     6. Generates post-processing plots
#
# Inputs:
#   None (all parameters are loaded from input_file.jl)
#
# Outputs:
#   solution     - Dict : Final solution state containing flow fields,
#                          mesh data, and solver metadata
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
using Arpack

## Solver root directory (for data file paths)
solver_dir = joinpath(@__DIR__, "..", "..") * "/"
solution = Dict{String, Any}()
solution["solver_dir"] = solver_dir

## Configure plot output directory (all figures saved here)
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output"))

## Load input configuration
filename = joinpath(@__DIR__, "input_file.jl")
LOAD_INPUT_VARIABLES(filename)

## Initialize chemistry
chemistry = SET_CHEMISTRY(solution)

## Initialize mesh, initial condition, and shock position
solution = INITIALIZATION(solution, solution_save, chemistry)

## Plot initial setup
PLOT_INITIAL_SET_UP(solution)

## Visualize velocity
plot_x_coords = solution["mesh"]["x"]
plot_y_coords = solution["mesh"]["y"]
plot_c_data   = solution["var"]["rho_u"][2:end-1, 2:end-1] ./ solution["var"]["rho"][2:end-1, 2:end-1]
(fig_vel, ax_vel) = curvilinear_heatmap(plot_x_coords, plot_y_coords, plot_c_data,
                    colorbar_label=L"$u$",
                    xlabel_str="x/h", ylabel_str="y/h",
                    fig_size=auto_figure_size(plot_x_coords, plot_y_coords))
fig_vel.tight_layout()
save_or_display_figure(fig_vel, name="channel_velocity")

## Linearize A
solution["boundary_conditions"]["periodic"] = true
(solution, L) = LINEARIZE_L(solution, chemistry)

## Compute eigenvalues of stability analysis
n_modes = 4 # search 4 most unstable eigenmodes
(d, V) = eigs(complex(L); nev=n_modes, sigma=0.00 + 0.25im) # Search close to region of reference
for i in 1:n_modes
    @printf("eigenvalue ")
    @printf("%i", i)
    @printf(" = ")
    @printf("%f + %fi", real(d[i]), imag(d[i]))
    @printf("\n")
end

## Visualize eigenmodes global modal analysis
SET_PLOT_DEFAULTS(output_dir=joinpath(@__DIR__, "output", "modal"))
mode = 1 # Mode to visualize
T_plot = 0 # Time in which disturbance is plotted
freestream_disturbances = false
PLOT_MODES(freestream_disturbances, L, solution, chemistry, V[:, mode], T_plot) # Plot modes
