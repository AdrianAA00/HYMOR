# HYMOR: HYpersonic MOdal/non-modal, and Receptivity

An open-source MATLAB/Julia framework for global modal, non-modal, and receptivity analysis of high-enthalpy hypersonic flows.

## Overview

HYMOR solves the compressible Navier--Stokes equations on structured curvilinear grids and provides integrated tools for linear stability analysis of high-enthalpy flows. The solver uses a shock-fitting formulation that treats the bow shock as a sharp discontinuity, reproducing the exact linear interaction analysis (LIA) response to infinitesimal disturbances and eliminating spurious artifacts associated with shock-capturing methods.

## Capabilities

### Flow Solver

- **Governing equations**: 2D and 3D-axisymmetric compressible Navier--Stokes equations discretized with a second-order finite-volume scheme.
- **Time integration**: Explicit fourth-order Runge--Kutta (RK4) and implicit backward Euler with dual time stepping.
- **Shock fitting**: Freely moving cubic-spline parametrization of the bow shock with Rankine--Hugoniot jump conditions, providing sub-grid accuracy and eliminating carbuncle-type instabilities.
- **Supported geometries**: Blunt wedge / blunt cone, cylinder / sphere, and wedge / cone. Custom curvilinear mappings can be added.
- **Boundary conditions**: Shock-compatibility (Rankine--Hugoniot), no-slip (adiabatic and isothermal), slip, Navier-slip, symmetry, non-reflecting characteristic (NRCBC), supersonic/subsonic inflow and outflow.

### Thermochemical Models

Five thermochemical models with increasing physical fidelity:

| Model | Description |
|-------|-------------|
| CPG | Calorically perfect gas |
| Frozen-RTV | Frozen chemistry with translational-rotational-vibrational equilibrium |
| Chemical-RTV | Chemical and translational-rotational-vibrational equilibrium |
| Chemical-RTVE | Chemical and translational-rotational-vibrational-electronic equilibrium |
| NonEq-RTVE | Chemical non-equilibrium with translational-rotational-vibrational-electronic equilibrium |

Equilibrium properties are computed via Cantera and stored as radial-basis-function fits for fast evaluation (~1000x speedup over online Cantera calls). Supported atmospheric compositions include Earth, Mars, and pure CO2.

Transport properties are available through Sutherland's law or collision-integral models (Wilke's mixing rule for viscosity, Mathur--Saxena rule for thermal conductivity).

### Linear Stability Analysis

**Global modal analysis** -- Solves the eigenvalue problem of the linearized Navier--Stokes operator to identify exponentially growing modes. Uses a GPU-accelerated exponential spectral transformation within an Implicitly Restarted Arnoldi Method (ARPACK), avoiding the dense LU/QR factorizations of the shift-and-invert approach and enabling analysis on large grids with minimal memory.

**Transient growth (non-modal) analysis** -- Computes the optimal initial disturbance that maximizes energy amplification over a prescribed time horizon, measured using Chu's energy norm. The optimization is formulated as a generalized Rayleigh quotient and solved via the Lanczos iteration with GPU-accelerated matrix exponentials.

**Freestream receptivity analysis** -- Determines the optimal freestream disturbance that maximizes post-shock energy amplification. Freestream perturbations are represented via a Fourier decomposition in the streamwise direction and coupled to the post-shock domain through the linearized Rankine--Hugoniot conditions. The solver automatically computes the optimal shape functions in the wall-normal direction for a user-specified set of streamwise wavenumbers.

All three analysis modes support automatic linearization of the discrete nonlinear operators (including shock perturbations in the shock-fitting formulation) and are compatible with all available thermochemical models.

## Repository Structure

```
HYMOR/
├── chemistry/                        Thermochemical models and property fits
│   ├── equilibrium_models/           Precomputed equilibrium property tables
│   ├── non_equilibrium_models/       Non-equilibrium relaxation parameters
│   └── shock_jump_properties/        Rankine--Hugoniot property data
├── tutorials/                        Example cases (see Tutorials section below)
│   ├── channel_flow/                 Incompressible channel flow (validation)
│   ├── hypersonic_cylinder/          Hypersonic flow over a cylinder
│   └── hypersonic_blunt_cone/        Entry capsule with real-gas effects
└── utils/                            Core solver modules
    ├── Energy_budgets/               Chu energy norm and budget decomposition
    ├── Initialization/               Solution setup and restart handling
    ├── Mesh/                         Curvilinear mesh generation
    ├── Operators/                    Flux computation and PDE discretization
    ├── Postprocessing/               Visualization routines
    ├── Shock_fitting/                Shock evolution and jump conditions
    ├── Stability_analysis/           Linear stability tools
    │   ├── Eigenvalues/              Arnoldi/Lanczos eigensolvers (GPU)
    │   ├── Modal_stability_analysis/ Linearization and modal analysis
    │   ├── Transient_growth_downstream/  Non-modal transient growth
    │   └── Freestream_receptivity/   Freestream receptivity optimization
    └── Time_marching/                Explicit and implicit time integrators
```

## Getting Started

### MATLAB Setup

#### Requirements

- **MATLAB R2023b** or later

#### Required Toolboxes

| Toolbox | Purpose |
|---------|---------|
| [Optimization Toolbox](https://www.mathworks.com/products/optimization.html) | `fsolve` and `optimoptions` for solving nonlinear Rankine-Hugoniot shock jump conditions |
| [Curve Fitting Toolbox](https://www.mathworks.com/products/curvefitting.html) | `csaps` and `ppval` for cubic smoothing spline fitting during shock resampling |

#### Optional Toolboxes

| Toolbox | Purpose |
|---------|---------|
| [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) | GPU acceleration via `gpuArray`, `gpuDevice`, and `gather` for eigenvalue solvers and stability analysis. CPU fallback is available if no GPU is present |

#### Optional External Software

| Software | Purpose |
|----------|---------|
| [Cantera](https://cantera.org/) | Only needed for generating new thermochemical property fits. Pre-computed property tables are included in `chemistry/equilibrium_models/` |

#### Verifying Your Installation

You can check that the required toolboxes are installed by running in the MATLAB Command Window:

```matlab
v = ver;
required = {'Optimization Toolbox', 'Curve Fitting Toolbox'};
for i = 1:length(required)
    if any(strcmp({v.Name}, required{i}))
        fprintf('%s: installed\n', required{i});
    else
        fprintf('%s: NOT FOUND\n', required{i});
    end
end
```

#### Quick Start

1. Open MATLAB and add the project to your path:
   ```matlab
   addpath(genpath('path/to/HYMOR'));
   ```
2. Navigate to a tutorial directory and run the driver script:
   ```matlab
   cd tutorials/channel_flow/
   run('main.m')
   ```

---

### Julia Setup

#### Requirements

- **Julia 1.9** or later

#### Installation

1. Navigate to the project root directory and start Julia:
   ```bash
   cd path/to/HYMOR
   julia
   ```

2. Activate the project environment and install all dependencies:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

   This will automatically download and install all required packages declared in `Project.toml`.

#### Dependencies

The Julia implementation depends on the following packages (resolved automatically by `Pkg.instantiate()`):

| Package | Purpose |
|---------|---------|
| Arpack | Eigenvalue computations via Implicitly Restarted Arnoldi Method |
| BSplineKit | B-spline curve fitting for shock resampling |
| CUDA | GPU acceleration and CUSPARSE sparse matrix operations |
| Interpolations | Numerical interpolation methods for thermochemical properties |
| JLD2 | HDF5-based data serialization for restart files |
| LinearMaps | Matrix-free linear operators for stability analysis |
| LoopVectorization | Performance optimization for computationally intensive loops |
| NaturalNeighbours | Scattered data interpolation for chemistry models |
| NLsolve | Nonlinear equation solving for Rankine-Hugoniot conditions |
| Plots | Visualization and postprocessing |
| PyPlot | Matplotlib-based plotting backend |
| LaTeXStrings | LaTeX-formatted labels in figures |

Standard library packages (`LinearAlgebra`, `SparseArrays`, `Printf`, `Dates`, `Statistics`) are also used and ship with Julia.

#### Quick Start

```julia
using Pkg
Pkg.activate(".")
include("tutorials/channel_flow/main.jl")
```

---

### Tutorials

Each tutorial directory contains driver scripts (`main.m` / `main.jl`), configuration files (`input_file.m` / `input_file.jl`), and a detailed README walkthrough. Click a tutorial below to get started:

| Tutorial | Description | Link |
|----------|-------------|------|
| Channel Flow | Incompressible channel flow stability analysis -- simplest starting point for validating the solver against reference eigenvalues and eigenmodes | [tutorials/channel_flow](tutorials/channel_flow/) |
| Hypersonic Cylinder | Hypersonic flow over a cylinder with shock fitting and modal stability analysis | [tutorials/hypersonic_cylinder](tutorials/hypersonic_cylinder/) |
| Hypersonic Entry Capsule | Entry capsule with real-gas thermochemistry, shock fitting, and full stability pipeline | [tutorials/hypersonic_blunt_cone](tutorials/hypersonic_blunt_cone/) |

To run any tutorial:

**MATLAB**
```matlab
cd tutorials/hypersonic_cylinder/
run('main.m')
```

**Julia**
```julia
using Pkg; Pkg.activate(".")
include("tutorials/hypersonic_cylinder/main.jl")
```

The driver script executes the full pipeline: initialization, base-flow computation, and (where included) modal stability analysis, transient growth, and freestream receptivity analysis.

### Configuring a Simulation

All simulation parameters are defined in the input file (`input_file.m` for MATLAB, `input_file.jl` for Julia). Key settings include:

- **Geometry and mesh**: `solution.curvilinear_mapping`, `solution.mesh`
- **Freestream conditions**: `solution.freestream` (velocity, density, temperature, Reynolds number)
- **Thermochemistry**: `solution.chemistry` (model type, atmospheric composition)
- **Shock fitting**: `solution.shock` (enable/disable, spline method, relaxation)
- **Time integration**: `solution.time_integration` (scheme, CFL, number of iterations)
- **Stability analysis**: `solution.stability_analysis` (eigensolver, perturbation magnitude)

## References

A. Anton-Alvarez and A. Lozano-Duran, "HYMOR: An open-source package for modal, non-modal,
  and receptivity analysis in high-enthalpy hypersonic vehicles."

## License

This project is released under the [MIT License](LICENSE).
