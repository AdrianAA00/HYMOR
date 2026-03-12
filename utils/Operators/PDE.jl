function PDE(s::Dict{String, Any}, chemistry::Dict{String, Any})
# PDE - Assemble all PDE flux contributions for a single time-advancement step.
#
# Orchestrates the full right-hand-side evaluation: updates equilibrium
# chemistry and thermodynamic properties, applies boundary conditions
# (including shock tracking), and computes inviscid and viscous flux
# contributions. Supports both 2D planar and 3D axisymmetric formulations.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract string comparisons once
    is_linearize = s["linearize"]::Bool
    shock_dict   = s["shock"]::Dict{String, Any}
    pde_dim      = s["PDE_dimension"]::String

    ## Update thermodynamic properties for new boundary conditions
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry)
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry)
    s = UPDATE_SOUND_SPEED(s, chemistry)

    ## Update shock position based on pressure signal from flow
    if shock_dict["enabled"]
        s = UPDATE_SHOCK_BC(s, chemistry)
    end

    ## Apply remaining boundary conditions
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry)

    ## Compute inviscid and viscous fluxes
    if pde_dim == "2D"
        s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry)
        s = FLUX_MASS(s)
        s = FLUX_MOMENTUM(s)
        s = FLUX_TOTAL_ENERGY(s)
        s = VISCOUS_FLUX(s)
    elseif pde_dim == "3D-axisymmetric"
        s = FLUX_NON_EQUILIBRIUM_CHEMISTRY(s, chemistry)
        s = FLUX_MASS_3D(s)
        s = FLUX_MOMENTUM_3D(s)
        s = FLUX_TOTAL_ENERGY_3D(s)
        s = VISCOUS_FLUX_3D(s)
    end

    return s
end
