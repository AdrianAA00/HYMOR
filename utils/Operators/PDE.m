function s = PDE(s, chemistry)
% PDE  Assemble all PDE flux contributions for a single time-advancement step.
%
%   s = PDE(s, chemistry)
%
%   Orchestrates the full right-hand-side evaluation: updates equilibrium
%   chemistry and thermodynamic properties, applies boundary conditions
%   (including shock tracking), and computes inviscid and viscous flux
%   contributions. Supports both 2D planar and 3D axisymmetric formulations.
%
%   Inputs:
%       s  - Solution struct containing the full flow state, mesh
%                   geometry, and configuration flags:
%                   .linearize  - If true, skip chemistry equilibrium update
%                   .PDE_dimension  - "2D" or "3D-axisymmetric"
%       chemistry - Chemistry model struct (lookup tables / evaluation fns)
%
%   Outputs:
%       s  - Updated s struct with all flux fields populated:
%                   .flux_rho, .flux_rho_u, .flux_rho_v, .flux_rho_E
%
%   Notes:
%       - In linearized mode, chemistry equilibrium is not recomputed.
%       - The 2D branch calls FLUX_MASS, FLUX_MOMENTUM, FLUX_TOTAL_ENERGY,
%         and VISCOUS_FLUX.
%       - The 3D-axisymmetric branch calls the _3D variants of each.
%
% Part of: Hypersonics Stability MATLAB Solver - Operators Module

    %% Update thermodynamic properties for new boundary conditions
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
    s = UPDATE_SOUND_SPEED(s, chemistry);

    %% Update shock position based on pressure signal from flow
    if s.shock.enabled
        s = UPDATE_SHOCK_BC(s, chemistry);
    end

    %% Apply remaining boundary conditions
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry);

    %% Compute inviscid and viscous fluxes
    if s.PDE_dimension == "2D"
        s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry);
        s = FLUX_MASS(s);
        s = FLUX_MOMENTUM(s);
        s = FLUX_TOTAL_ENERGY(s);
        s = VISCOUS_FLUX(s);
    elseif s.PDE_dimension == "3D-axisymmetric"
        s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry);
        s = FLUX_MASS_3D(s);
        s = FLUX_MOMENTUM_3D(s);
        s = FLUX_TOTAL_ENERGY_3D(s);
        s = VISCOUS_FLUX_3D(s);
    end
end
