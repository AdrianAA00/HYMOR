function UPDATE_PRESSURE(rho, T, chemistry::Dict{String, Any}, s::Dict{String, Any})
# UPDATE_PRESSURE - Compute pressure from density and temperature using chemistry tables.
#
# Evaluates the pressure using the equation of state p = rho * e * (gamma* - 1).
#
# Note: The original MATLAB code references s.chemistry.is_chemistry_enabled
# from the caller's workspace. In Julia, we pass s explicitly.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    chem = s["chemistry"]::Dict{String, Any}

    ## Compute pressure from equation of state
    if chem["is_chemistry_enabled"]::Bool
        e = chemistry["eval_e"](T, rho)
        gamma_star = chemistry["eval_gamma_star"](rho, e)
        p = @. rho * e * (gamma_star - 1)
    else
        p = zeros(Float64, size(rho))  # Placeholder for non-chemistry case
    end

    return p
end
