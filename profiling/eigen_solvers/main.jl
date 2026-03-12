# main.jl - Entry point for eigensolver scaling profiling
#
# Description:
#   Test scaling and performance of eigenvalue, transient growth and
#   freestream receptivity solvers across multiple grid sizes.
#   Timings are printed to screen and saved to the "timings/" folder.
#
# Inputs:
#   None (all parameters are loaded from input_file.jl)
#
# Outputs:
#   Timing text files in the "timings/" sub-directory.
#
# Usage:
#   Run this script from the test case directory:
#       julia main.jl
#
# Part of: Hypersonics Stability MATLAB Solver - Test Cases Module

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

## Create timings output directory
timings_dir = joinpath(@__DIR__, "timings")
mkpath(timings_dir)

## Load input configuration
filename = joinpath(@__DIR__, "input_file.jl")
LOAD_INPUT_VARIABLES(filename)

## Load chemistry model
if !solution["restart"]
    chemistry = SET_CHEMISTRY(solution)
end

## Run simulation (256 x 64 grid) for base flow
solution["restart"]                        = false
solution["time_integration"]["N_iter"]     = 100
solution["mesh"]["Nchi"]                   = 2^8
solution["mesh"]["Neta"]                   = 2^6
solution["freestream"]["u"]                = 2000       # [m/s]
solution["freestream"]["Re"]               = 2000       # Reynolds number
solution_save = nothing
solution = INITIALIZATION(solution, solution_save, chemistry)

solution_base = deepcopy(solution)
t_start = time()
solution = RUN_SIMULATION(solution, solution_base, chemistry)
t_elapsed = time() - t_start
println("Base-flow elapsed time: $t_elapsed seconds")
solution_save = deepcopy(solution)

# Increase Reynolds for tests
solution["freestream"]["Re"] = 50000

## ========================================================================
#  Grid-scaling timing study
#  ========================================================================

# Define grid sizes to test (Nchi x Neta)
grid_Nchi = [2^8,  2^9,  2^10, 2^11, 2^12]
grid_Neta = [2^6,  2^7,  2^8,  2^9,  2^10]
n_grids   = 3; ## Number of grid sizes to test (use first n_grids entries from grid_Nchi and grid_Neta)
n_modes   = 10
N_l       = 40

## Warm-up pass (smallest grid) — triggers JIT compilation so timed runs are fair
println("Warm-up: compiling solver functions on smallest grid...")
solution["mesh"]["Nchi"] = grid_Nchi[1]
solution["mesh"]["Neta"] = grid_Neta[1]
solution["restart"]      = true
solution = INITIALIZATION(solution, solution_save, chemistry)
(solution, L_wu) = LINEARIZE_L(solution, chemistry)
(_, _) = EIGENVALUES(L_wu, solution, n_modes)
(_, _, _) = TRANSIENT_GROWTH_DOWNSTREAM(L_wu, solution, n_modes, 1.0)
w_infty_wu = 2*pi*collect(range(0, N_l, length=N_l+1))
(solution, L__wu) = LINEARIZE_L_(L_wu, solution, chemistry, w_infty_wu)
(_, _, _) = FREESTREAM_RECEPTIVITY(L__wu, solution, n_modes, 1.0, w_infty_wu)
L_wu = nothing; L__wu = nothing; GC.gc()
println("Warm-up complete.\n")

# Pre-allocate timing arrays
t_eigenvalues  = zeros(n_grids)
t_tgd          = zeros(n_grids)
t_freestream   = zeros(n_grids)

# Print header
println()
println("==========================================================================")
println("  Eigensolver scaling study — Julia")
println("==========================================================================")
@printf("%-6s  %-6s  %-12s  %-16s  %-16s  %-16s\n",
        "Nchi", "Neta", "Nchi*Neta", "EIGENVALUES [s]", "TGD [s]", "FREESTREAM [s]")
println("--------------------------------------------------------------------------")

# Write file headers only if files do not exist yet
eig_file = joinpath(timings_dir, "eigenvalues_timings_julia.txt")
tgd_file = joinpath(timings_dir, "transient_growth_downstream_timings_julia.txt")
fr_file  = joinpath(timings_dir, "freestream_receptivity_timings_julia.txt")

if !isfile(eig_file)
    open(eig_file, "w") do fid
        write(fid, "# Eigensolver scaling: EIGENVALUES\n")
        @printf(fid, "# %-10s  %-10s  %-14s  %-10s  %-16s\n", "Nchi", "Neta", "Nchi*Neta", "Language", "Time [s]")
    end
end
if !isfile(tgd_file)
    open(tgd_file, "w") do fid
        write(fid, "# Eigensolver scaling: TRANSIENT_GROWTH_DOWNSTREAM\n")
        @printf(fid, "# %-10s  %-10s  %-14s  %-10s  %-16s\n", "Nchi", "Neta", "Nchi*Neta", "Language", "Time [s]")
    end
end
if !isfile(fr_file)
    open(fr_file, "w") do fid
        write(fid, "# Eigensolver scaling: FREESTREAM_RECEPTIVITY\n")
        @printf(fid, "# %-10s  %-10s  %-14s  %-10s  %-16s\n", "Nchi", "Neta", "Nchi*Neta", "Language", "Time [s]")
    end
end

for ig in 1:n_grids
    global solution, t_start

    Nchi_i = grid_Nchi[ig]
    Neta_i = grid_Neta[ig]

    ## Set up grid and linearize
    solution["mesh"]["Nchi"] = Nchi_i
    solution["mesh"]["Neta"] = Neta_i
    solution["restart"]      = true
    solution = INITIALIZATION(solution, solution_save, chemistry)
    (solution, L) = LINEARIZE_L(solution, chemistry)

    ## Time EIGENVALUES
    t_start = time()
    (_, _) = EIGENVALUES(L, solution, n_modes)
    t_eigenvalues[ig] = time() - t_start

    ## Time TRANSIENT_GROWTH_DOWNSTREAM
    T_TGD = 1.0
    t_start = time()
    (_, _, _) = TRANSIENT_GROWTH_DOWNSTREAM(L, solution, n_modes, T_TGD)
    t_tgd[ig] = time() - t_start

    ## Linearize extended system
    w_infty = 2*pi*collect(range(0, N_l, length=N_l+1))
    (solution, L_) = LINEARIZE_L_(L, solution, chemistry, w_infty)

    ## Time FREESTREAM_RECEPTIVITY
    T_TGF = 1.0
    t_start = time()
    (_, _, _) = FREESTREAM_RECEPTIVITY(L_, solution, n_modes, T_TGF, w_infty)
    t_freestream[ig] = time() - t_start

    ## Print row
    @printf("%-6d  %-6d  %-12d  %-16.4f  %-16.4f  %-16.4f\n",
            Nchi_i, Neta_i, Nchi_i*Neta_i,
            t_eigenvalues[ig], t_tgd[ig], t_freestream[ig])

    ## Append this iteration's results to timing files
    open(eig_file, "a") do fid
        @printf(fid, "  %-10d  %-10d  %-14d  %-10s  %-16.6f\n",
                Nchi_i, Neta_i, Nchi_i*Neta_i, "Julia", t_eigenvalues[ig])
    end
    open(tgd_file, "a") do fid
        @printf(fid, "  %-10d  %-10d  %-14d  %-10s  %-16.6f\n",
                Nchi_i, Neta_i, Nchi_i*Neta_i, "Julia", t_tgd[ig])
    end
    open(fr_file, "a") do fid
        @printf(fid, "  %-10d  %-10d  %-14d  %-10s  %-16.6f\n",
                Nchi_i, Neta_i, Nchi_i*Neta_i, "Julia", t_freestream[ig])
    end

    L = nothing; L_ = nothing;
end

println("==========================================================================")
println()
println("Timing files saved in: ", timings_dir)
