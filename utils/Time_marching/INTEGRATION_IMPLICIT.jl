# INTEGRATION_IMPLICIT  Perform one implicit Euler update with relaxation.
#
#   (solution_temp, residual) = INTEGRATION_IMPLICIT(s, solution_temp)
#
#   Computes a forward-Euler predictor for each conserved variable (rho,
#   rho_u, rho_v, rho_E), evaluates the maximum absolute difference
#   (residual) between the current iterate and the predictor, and then
#   blends the two using the adaptive relaxation factor stored in
#   s["relax_factor_variable"]. The shock position is also updated
#   with an explicit Euler step.
#
#   Inputs:
#       s              (Dict{String,Any}) - Reference (time-level n) flow state
#                                           containing conservative variables
#                                           (s["var"].*), dt (s["time_integration"]["dt"]),
#                                           flow_cells mask (s["shock"]["flow_cells"]),
#                                           shock data (s["shock"].*), and the current
#                                           relaxation factor (s["relax_factor_variable"]).
#       solution_temp  (Dict{String,Any}) - Current implicit iterate whose fluxes
#                                           have already been evaluated by PDE.
#
#   Outputs:
#       solution_temp  (Dict{String,Any}) - Updated iterate after relaxation blending
#                                           for rho, rho_u, rho_v, rho_E and after
#                                           an explicit shock-position update.
#       residual       (Float64)          - Maximum absolute residual across all four
#                                           conserved variables (excluding shock and
#                                           ghost cells).
#
#   Notes:
#       - For shock-fitted runs the fluxes are masked by
#         s["shock"]["flow_cells"] and the residual is zeroed above the shock
#         row to avoid polluting the convergence measure.
#       - The relaxation update is:
#           q^{k+1} = relax * q^{k} + (1 - relax) * q_predictor
#         where relax = s["relax_factor_variable"].
#       - Shock position is advanced with an explicit Euler step scaled
#         by s["shock"]["relaxation"].
#
# Part of: Hypersonics Stability Julia Solver - Time Marching Module

function INTEGRATION_IMPLICIT(s::Dict{String,Any}, solution_temp::Dict{String,Any})

    ## Extract dict refs with type assertions for compiler optimization
    s_var       = s["var"]::Dict{String,Any}
    s_shock     = s["shock"]::Dict{String,Any}
    s_ti        = s["time_integration"]::Dict{String,Any}
    st_var      = solution_temp["var"]::Dict{String,Any}
    st_flux     = solution_temp["flux"]::Dict{String,Any}
    st_shock    = solution_temp["shock"]::Dict{String,Any}

    s_rho       = s_var["rho"]::Matrix{Float64}
    s_rho_u     = s_var["rho_u"]::Matrix{Float64}
    s_rho_v     = s_var["rho_v"]::Matrix{Float64}
    s_rho_E     = s_var["rho_E"]::Matrix{Float64}

    st_rho      = st_var["rho"]::Matrix{Float64}
    st_rho_u    = st_var["rho_u"]::Matrix{Float64}
    st_rho_v    = st_var["rho_v"]::Matrix{Float64}
    st_rho_E    = st_var["rho_E"]::Matrix{Float64}

    fl_rho      = st_flux["rho"]::Matrix{Float64}
    fl_rho_u    = st_flux["rho_u"]::Matrix{Float64}
    fl_rho_v    = st_flux["rho_v"]::Matrix{Float64}
    fl_rho_E    = st_flux["rho_E"]::Matrix{Float64}

    dt           = s_ti["dt"]::Float64
    shock_en     = s_shock["enabled"]::Bool
    relax        = s["relax_factor_variable"]::Float64
    one_m_relax  = 1.0 - relax

    flow_cells   = s_shock["flow_cells"]::Matrix{Float64}
    cell_indices = s_shock["cell_indices"]::Matrix{Int}
    Nchi         = s["mesh"]["Nchi"]::Int

    @views begin

    ## Density (rho) update
    if shock_en
        temp = @. s_rho[2:end-1, 2:end-1] + dt * fl_rho * flow_cells
    else
        temp = @. s_rho[2:end-1, 2:end-1] + dt * fl_rho
    end
    residual_temp = @. abs(st_rho[2:end-1, 3:end-2] - temp[:, 2:end-1])
    if shock_en
        @inbounds for i in 1:Nchi
            residual_temp[i, cell_indices[i, 1]-1:end] .= 0.0
        end
    end
    res_rho = maximum(residual_temp)
    @. st_rho[2:end-1, 2:end-1] = st_rho[2:end-1, 2:end-1] * relax + one_m_relax * temp

    ## X-momentum (rho_u) update
    if shock_en
        temp = @. s_rho_u[2:end-1, 2:end-1] + dt * fl_rho_u * flow_cells
    else
        temp = @. s_rho_u[2:end-1, 2:end-1] + dt * fl_rho_u
    end
    residual_temp = @. abs(st_rho_u[2:end-1, 3:end-2] - temp[:, 2:end-1])
    if shock_en
        @inbounds for i in 1:Nchi
            residual_temp[i, cell_indices[i, 1]-1:end] .= 0.0
        end
    end
    res_rho_u = maximum(residual_temp)
    @. st_rho_u[2:end-1, 2:end-1] = st_rho_u[2:end-1, 2:end-1] * relax + one_m_relax * temp

    ## Y-momentum (rho_v) update
    if shock_en
        temp = @. s_rho_v[2:end-1, 2:end-1] + dt * fl_rho_v * flow_cells
    else
        temp = @. s_rho_v[2:end-1, 2:end-1] + dt * fl_rho_v
    end
    residual_temp = @. abs(st_rho_v[2:end-1, 3:end-2] - temp[:, 2:end-1])
    if shock_en
        @inbounds for i in 1:Nchi
            residual_temp[i, cell_indices[i, 1]-1:end] .= 0.0
        end
    end
    res_rho_v = maximum(residual_temp)
    @. st_rho_v[2:end-1, 2:end-1] = st_rho_v[2:end-1, 2:end-1] * relax + one_m_relax * temp

    ## Total energy (rho_E) update
    if shock_en
        temp = @. s_rho_E[2:end-1, 2:end-1] + dt * fl_rho_E * flow_cells
    else
        temp = @. s_rho_E[2:end-1, 2:end-1] + dt * fl_rho_E
    end
    residual_temp = @. abs(st_rho_E[2:end-1, 3:end-2] - temp[:, 2:end-1])
    if shock_en
        @inbounds for i in 1:Nchi
            residual_temp[i, cell_indices[i, 1]-1:end] .= 0.0
        end
    end
    res_rho_E = maximum(residual_temp)
    @. st_rho_E[2:end-1, 2:end-1] = st_rho_E[2:end-1, 2:end-1] * relax + one_m_relax * temp

    ## Overall residual (avoid allocating a Vector)
    residual = max(res_rho, res_rho_u, res_rho_v, res_rho_E)

    ## Update shock position (explicit Euler step)
    s_pts_x   = s_shock["points_x"]::Matrix{Float64}
    s_pts_y   = s_shock["points_y"]::Matrix{Float64}
    st_spd_x  = st_shock["speed_x"]::Matrix{Float64}
    st_spd_y  = st_shock["speed_y"]::Matrix{Float64}
    shock_relax = s_shock["relaxation"]::Float64
    st_shock["points_x"] = @. s_pts_x + st_spd_x * dt * shock_relax
    st_shock["points_y"] = @. s_pts_y + st_spd_y * dt * shock_relax

    end # @views

    return solution_temp, residual
end
