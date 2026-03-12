function MIN_RHO(s::Dict{String, Any})
# MIN_RHO - Enforce minimum density threshold on the solution field.
#
# Clamps the density field to a prescribed minimum value within the
# interior cells. If any interior density values fall below the minimum
# threshold, they are replaced and a warning message is printed.
#
# Part of: Hypersonics Stability Julia Solver - Operators Module

    ## Extract dict refs
    var     = s["var"]::Dict{String, Any}
    rho     = var["rho"]::Matrix{Float64}
    rho_min = s["numerical_dissipation"]["rho_min"]

    ## Apply minimum density threshold only to interior
    @views interior = rho[2:end-1, 2:end-1]
    limited = false
    @inbounds for j in axes(interior, 2), i in axes(interior, 1)
        if interior[i, j] < rho_min
            interior[i, j] = rho_min
            limited = true
        end
    end

    if limited
        println()
        println(" Rho limited activated ")
        println()
    end

    return s
end
