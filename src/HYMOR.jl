# HYMOR.jl - Main module for HYMOR - HYpersonic MOdal/non-modal, and Receptivity
#
# Wraps all solver utilities, chemistry models, mesh generation, operators,
# shock fitting, time marching, stability analysis, and postprocessing into
# a single Julia module with proper dependency ordering.
#
# Usage:
#   import Pkg; Pkg.activate("path/to/Hypersonics_stability_MATLAB")
#   using HYMOR

module HYMOR

# ── External dependencies ─────────────────────────────────────────────
using LinearAlgebra
using SparseArrays
using Printf
using Dates
using Statistics

using Arpack
using BSplineKit
using CUDA
using Interpolations
using LinearMaps
using JLD2
using LaTeXStrings
using LoopVectorization
using NaturalNeighbours
using NLsolve
using Plots
import PyPlot

# ── Exports ───────────────────────────────────────────────────────────
# Chemistry
export ATMOSPHERE_PROPERTIES, EARTH_ATMOSPHERE, MARS_ATMOSPHERE
export READ_GAS_PROPERTIES
export READ_NON_EQUILIBRIUM_MODEL
export READ_SHOCK_JUMP_CONDITIONS
export SOLVE_RANKINE_HUGONIOT_CHEMISTRY
export SET_CHEMISTRY
export TEST_RANKINE_HUGONIOT_SOLVER, PLOT_TEST_RANKINE_HUGONIOT

# Mesh
export GO_TO_ELEMENT_SPACE, GO_TO_PHYSICAL_SPACE
export EXTEND_TO_GHOST_POINTS, GET_EXTENDED_CELL_CENTERS
export GET_CELL_CORNERS_SMOOTH, GENERATE_MESH
export CREATE_MASKS, CHECK_MESH

# Operators – core
export DERIVATIVE, DERIVATIVE_EXT
export FILTERING, NUMERICAL_VISCOSITY
export MIN_RHO, REFINEMENT
# Operators – thermodynamics & BCs
export UPDATE_PRESSURE, UPDATE_SOUND_SPEED
export UPDATE_CHEMISTRY_EQUILIBRIUM, UPDATE_THERMODYNAMIC_PROPERTIES
export APPLY_BOUNDARY_CONDITIONS
# Operators – fluxes
export FLUX_MASS, FLUX_MOMENTUM, FLUX_TOTAL_ENERGY
export FLUX_NON_EQUILIBRIUM_CHEMISTRY, VISCOUS_FLUX
export FLUX_MASS_3D, FLUX_MOMENTUM_3D, FLUX_TOTAL_ENERGY_3D, VISCOUS_FLUX_3D
# Operators – PDE
export NON_LINEAR_DYNAMICS_NO_DISCONTINUITY, PDE

# Shock fitting
export _spline_derivative
export COMPUTE_UPSTREAM_CONDITIONS, COMPUTE_BETA
export DETECT_CELLS_SHOCKED, EXTRAPOLATE_CELLS_SHOCK
export LEAST_SQUARES_SHOCK_POINTS
export UPDATE_FIELD_UPSTREAM, UPDATE_FLOW_CELLS
export INITIAL_SOLUTION_POST_SHOCK, INITIAL_SOLUTION_POST_SHOCK_MODIFIED
export MOVE_SHOCK, REMOVE_AUX_SHOCK_CELLS
export UPDATE_SHOCK_JUMP_PROPERTIES, UPDATE_SHOCK_BC
export CHECK_MESH_FOLLOWS_SHOCK

# Initialization
export LOAD_INPUT_VARIABLES
export VERIFY_INPUT_DATA, VARIABLES_INITIALIZATION
export SET_FREESTREAM_PROPERTIES, SHOCK_LINE_INITIALIZATION
export INITIAL_SOLUTION
export RESTART_SOLUTION, RESTART_FROM_FILE, START_SOLUTION
export INITIALIZATION

# Time marching
export CFL_TIMESTEP, INTEGRATION
export INTEGRATION_IMPLICIT
export ADVECT_FLOW_EXPLICIT, ADVECT_FLOW_IMPLICIT, ADVECT_FLOW
export RUN_SIMULATION

# Energy budgets
export GET_ENERGY_NORM, GET_CHU_COMPONENTS
export COMPUTE_BUDGETS_KINETIC, COMPUTE_BUDGETS_ENTROPY, COMPUTE_BUDGETS_SHOCK
export ENERGY_INFLOW

# Postprocessing
export SET_PLOT_DEFAULTS, auto_figure_size, curvilinear_heatmap, save_or_display_figure
export PLOT_INITIAL_SET_UP, PLOT_INFLECTION_POINTS
export PLOT_RUNNING, PLOT_MODES, PLOT_MODES_CHU_NORM
export ppval

# Stability analysis – eigenvalues
export ARNOLDI_POWER_GPU, EIGENVALUES
# Stability analysis – modal
export LINEARIZE_EQUATIONS_NO_DISCONTINUITY, CHECK_LINEARIZATION_L, LINEARIZE_L
# Stability analysis – transient growth downstream
export CONSTRUCT_M, CONSTRUCT_M_HAT, CONSTRUCT_C, CONSTRUCT_R
export TRANSIENT_GROWTH_DOWNSTREAM
# Stability analysis – freestream receptivity
export GET_T_REF, CONSTRUCT_R_, CONSTRUCT_M_, CONSTRUCT_M_INFTY_
export LINEARIZE_L_, LINEARIZE_B, CHECK_LINEARIZATION_B
export LINEAR_INTEGRATION_AND_GAINS, FREESTREAM_RECEPTIVITY

# ── MATLAB compatibility helpers ────────────────────────────────────
"""Evaluate a BSplineKit spline (or any callable) at the given points — mirrors MATLAB's ppval."""
ppval(spline, x) = spline.(x)

# ── Layer 1: Chemistry utilities (no internal dependencies) ──────────
include("../chemistry/ATMOSPHERE_PROPERTIES.jl")
include("../chemistry/READ_GAS_PROPERTIES.jl")
include("../chemistry/READ_NON_EQUILIBRIUM_MODEL.jl")
include("../chemistry/READ_SHOCK_JUMP_CONDITIONS.jl")

# ── Layer 2: Chemistry solvers ───────────────────────────────────────
include("../chemistry/SOLVE_RANKINE_HUGONIOT_CHEMISTRY.jl")
include("../chemistry/SET_CHEMISTRY.jl")
include("../chemistry/TEST_RANKINE_HUGONIOT_SOLVER.jl")

# ── Layer 3: Mesh generation ─────────────────────────────────────────
include("../utils/Mesh/GO_TO_ELEMENT_SPACE.jl")
include("../utils/Mesh/GO_TO_PHYSICAL_SPACE.jl")
include("../utils/Mesh/EXTEND_TO_GHOST_POINTS.jl")
include("../utils/Mesh/GET_EXTENDED_CELL_CENTERS.jl")
include("../utils/Mesh/GET_CELL_CORNERS_SMOOTH.jl")
include("../utils/Mesh/GENERATE_MESH.jl")
include("../utils/Mesh/CREATE_MASKS.jl")
include("../utils/Mesh/CHECK_MESH.jl")

# ── Layer 4: Operators – core spatial utilities ──────────────────────
include("../utils/Operators/DERIVATIVE.jl")
include("../utils/Operators/DERIVATIVE_EXT.jl")
include("../utils/Operators/FILTERING.jl")
include("../utils/Operators/NUMERICAL_VISCOSITY.jl")
include("../utils/Operators/MIN_RHO.jl")
include("../utils/Operators/REFINEMENT.jl")

# ── Layer 5: Operators – thermodynamics ──────────────────────────────
include("../utils/Operators/UPDATE_PRESSURE.jl")
include("../utils/Operators/UPDATE_SOUND_SPEED.jl")
include("../utils/Operators/UPDATE_CHEMISTRY_EQUILIBRIUM.jl")
include("../utils/Operators/UPDATE_THERMODYNAMIC_PROPERTIES.jl")

# ── Layer 6: Operators – boundary conditions ─────────────────────────
include("../utils/Operators/APPLY_BOUNDARY_CONDITIONS.jl")

# ── Layer 7: Operators – flux computations ───────────────────────────
include("../utils/Operators/FLUX_MASS.jl")
include("../utils/Operators/FLUX_MOMENTUM.jl")
include("../utils/Operators/FLUX_TOTAL_ENERGY.jl")
include("../utils/Operators/FLUX_NON_EQUILIBRIUM_CHEMISTRY.jl")
include("../utils/Operators/VISCOUS_FLUX.jl")
include("../utils/Operators/FLUX_MASS_3D.jl")
include("../utils/Operators/FLUX_MOMENTUM_3D.jl")
include("../utils/Operators/FLUX_TOTAL_ENERGY_3D.jl")
include("../utils/Operators/VISCOUS_FLUX_3D.jl")

# ── Layer 8: Operators – PDE assembly ────────────────────────────────
include("../utils/Operators/PDE.jl")
include("../utils/Operators/NON_LINEAR_DYNAMICS_NO_DISCONTINUITY.jl")

# ── Layer 9: Shock fitting ───────────────────────────────────────────
include("../utils/Shock_fitting/CSAPS.jl")
include("../utils/Shock_fitting/COMPUTE_UPSTREAM_CONDITIONS.jl")
include("../utils/Shock_fitting/COMPUTE_BETA.jl")
include("../utils/Shock_fitting/DETECT_CELLS_SHOCKED.jl")
include("../utils/Shock_fitting/EXTRAPOLATE_CELLS_SHOCK.jl")
include("../utils/Shock_fitting/LEAST_SQUARES_SHOCK_POINTS.jl")
include("../utils/Shock_fitting/UPDATE_FIELD_UPSTREAM.jl")
include("../utils/Shock_fitting/UPDATE_FLOW_CELLS.jl")
include("../utils/Shock_fitting/INITIAL_SOLUTION_POST_SHOCK.jl")
include("../utils/Shock_fitting/INITIAL_SOLUTION_POST_SHOCK_MODIFIED.jl")
include("../utils/Shock_fitting/MOVE_SHOCK.jl")
include("../utils/Shock_fitting/REMOVE_AUX_SHOCK_CELLS.jl")
include("../utils/Shock_fitting/UPDATE_SHOCK_JUMP_PROPERTIES.jl")
include("../utils/Shock_fitting/UPDATE_SHOCK_BC.jl")
include("../utils/Shock_fitting/CHECK_MESH_FOLLOWS_SHOCK.jl")

# ── Layer 10: Initialization ─────────────────────────────────────────
include("../utils/Initialization/LOAD_INPUT_VARIABLES.jl")
include("../utils/Initialization/VERIFY_INPUT_DATA.jl")
include("../utils/Initialization/VARIABLES_INITIALIZATION.jl")
include("../utils/Initialization/SET_FREESTREAM_PROPERTIES.jl")
include("../utils/Initialization/SHOCK_LINE_INITIALIZATION.jl")
include("../utils/Initialization/INITIAL_SOLUTION.jl")
include("../utils/Initialization/RESTART_SOLUTION.jl")
include("../utils/Initialization/RESTART_FROM_FILE.jl")
include("../utils/Initialization/START_SOLUTION.jl")
include("../utils/Initialization/INITIALIZATION.jl")

# ── Layer 11: Time marching ──────────────────────────────────────────
include("../utils/Time_marching/CFL_TIMESTEP.jl")
include("../utils/Time_marching/INTEGRATION.jl")
include("../utils/Time_marching/INTEGRATION_IMPLICIT.jl")
include("../utils/Time_marching/ADVECT_FLOW_EXPLICIT.jl")
include("../utils/Time_marching/ADVECT_FLOW_IMPLICIT.jl")
include("../utils/Time_marching/ADVECT_FLOW.jl")
include("../utils/Time_marching/RUN_SIMULATION.jl")

# ── Layer 12: Energy budgets ─────────────────────────────────────────
include("../utils/Energy_budgets/GET_ENERGY_NORM.jl")
include("../utils/Energy_budgets/GET_CHU_COMPONENTS.jl")
include("../utils/Energy_budgets/COMPUTE_BUDGETS_KINETIC.jl")
include("../utils/Energy_budgets/COMPUTE_BUDGETS_ENTROPY.jl")
include("../utils/Energy_budgets/COMPUTE_BUDGETS_SHOCK.jl")
include("../utils/Energy_budgets/ENERGY_INFLOW.jl")

# ── Layer 13: Postprocessing ─────────────────────────────────────────
include("../utils/Postprocessing/PLOT_DEFAULTS.jl")
include("../utils/Postprocessing/PLOT_INITIAL_SET_UP.jl")
include("../utils/Postprocessing/PLOT_INFLECTION_POINTS.jl")
include("../utils/Postprocessing/PLOT_RUNNING.jl")
include("../utils/Postprocessing/PLOT_MODES.jl")
include("../utils/Postprocessing/PLOT_MODES_CHU_NORM.jl")

# ── Layer 14: Stability analysis ─────────────────────────────────────
# Eigenvalue solvers
include("../utils/Stability_analysis/Eigenvalues/ARNOLDI_POWER_GPU.jl")
include("../utils/Stability_analysis/Eigenvalues/EIGENVALUES.jl")

# Modal stability analysis
include("../utils/Stability_analysis/Modal_stability_analysis/LINEARIZE_EQUATIONS_NO_DISCONTINUITY.jl")
include("../utils/Stability_analysis/Modal_stability_analysis/CHECK_LINEARIZATION_L.jl")
include("../utils/Stability_analysis/Modal_stability_analysis/LINEARIZE_L.jl")

# Transient growth downstream
include("../utils/Stability_analysis/Transient_growth_downstream/CONSTRUCT_M.jl")
include("../utils/Stability_analysis/Transient_growth_downstream/CONSTRUCT_M_HAT.jl")
include("../utils/Stability_analysis/Transient_growth_downstream/CONSTRUCT_C.jl")
include("../utils/Stability_analysis/Transient_growth_downstream/CONSTRUCT_R.jl")
include("../utils/Stability_analysis/Transient_growth_downstream/TRANSIENT_GROWTH_DOWNSTREAM.jl")

# Freestream receptivity
include("../utils/Stability_analysis/Freestream_receptivity/GET_T_REF.jl")
include("../utils/Stability_analysis/Freestream_receptivity/CONSTRUCT_R_.jl")
include("../utils/Stability_analysis/Freestream_receptivity/CONSTRUCT_M_.jl")
include("../utils/Stability_analysis/Freestream_receptivity/CONSTRUCT_M_INFTY_.jl")
include("../utils/Stability_analysis/Freestream_receptivity/LINEARIZE_L_.jl")
include("../utils/Stability_analysis/Freestream_receptivity/LINEARIZE_B.jl")
include("../utils/Stability_analysis/Freestream_receptivity/CHECK_LINEARIZATION_B.jl")
include("../utils/Stability_analysis/Freestream_receptivity/LINEAR_INTEGRATION_AND_GAINS.jl")
include("../utils/Stability_analysis/Freestream_receptivity/FREESTREAM_RECEPTIVITY.jl")

end # module HYMOR
