# main.jl - Entry point for 3D axisymmetric blunted cone simulation
#
# Description:
#   Script to measure scaling of the time-marching solver without chemistry and no shock-fitting
#
# Usage:
#   Run this script from the test case directory:
#       julia main.jl
#
# Part of: Hypersonics Stability MATLAB Solver - Profiling Module

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

## Open output file for timing results
results_file = open(joinpath(@__DIR__, "timings", "timing_results_julia.txt"), "w")
println(results_file, "# Timing results - Julia")
println(results_file, "# Grid (Nchi x Neta)    Elapsed time (s)")
@printf(results_file, "%-25s %s\n", "Grid", "Elapsed_time_s")

## Time marching 100 x 100 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 100   # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 100   # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

# Warmup call (triggers compilation, result discarded)
solution_base = solution
RUN_SIMULATION(solution, solution_base, chemistry)
solution = INITIALIZATION(solution, solution_save, chemistry)
solution_base = nothing

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

## Time marching 200 x 200 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 200   # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 200   # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

## Time marching 400 x 400 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 400   # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 400   # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

## Time marching 800 x 800 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 800   # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 800   # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

## Time marching 1600 x 1600 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 1600  # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 1600  # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

## Time marching 3200 x 3200 grid, 1000 iterations
solution["restart"]                         = false
solution["time_integration"]["N_iter"]      = 1000  # Number of time steps
solution["mesh"]["Nchi"]                    = 3200  # Grid points along chi (streamwise)
solution["mesh"]["Neta"]                    = 3200  # Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = solution
elapsed = @elapsed solution = RUN_SIMULATION(solution, solution_base, chemistry)
@printf("Grid: %d x %d - Elapsed time: %.6f s\n", solution["mesh"]["Nchi"], solution["mesh"]["Neta"], elapsed)
@printf(results_file, "%-25s %.6f\n", "$(solution["mesh"]["Nchi"])x$(solution["mesh"]["Neta"])", elapsed)
solution_base = nothing

close(results_file)
