using NLsolve
using Printf

"""
    SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, w_1, rho_1, e_1; InitialGuess=[8.0, 1.2])

Vectorized Rankine-Hugoniot shock solver with real-gas chemistry.

Solves N independent 2x2 nonlinear systems simultaneously using a
single NLsolve call. The key optimisation is a single vectorised
chemistry evaluation per solver iteration instead of N separate calls.

# Arguments
- `chemistry::Dict{String,Any}`: Fitted chemistry structure from FIT_CHEMISTRY
  (must contain `"eval_gamma_star"`)
- `w_1`: Upstream velocity [m/s] (scalar or N-vector)
- `rho_1`: Upstream density [kg/m^3] (scalar or N-vector)
- `e_1`: Upstream specific internal energy [J/kg] (scalar or N-vector)

# Keyword Arguments
- `InitialGuess::Vector{Float64}`: [p0, q0] initial density and energy ratios (default: [8.0, 1.2])

# Returns
- `(rho_2, e_2, w_2)`: Post-shock density, energy, and velocity vectors

# Part of: Hypersonics Stability MATLAB Solver - Chemistry Module (Julia port)
"""
function SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry::Dict{String,Any},
                                          w_1, rho_1, e_1;
                                          InitialGuess::Vector{Float64}=[8.0, 1.2])

    ## Prepare input vectors
    w_1_vec   = vec(collect(Float64, w_1))
    rho_1_vec = vec(collect(Float64, rho_1))
    e_1_vec   = vec(collect(Float64, e_1))
    N = length(w_1_vec)

    if length(rho_1_vec) == 1
        rho_1_vec = fill(rho_1_vec[1], N)
    end
    if length(e_1_vec) == 1
        e_1_vec = fill(e_1_vec[1], N)
    end

    ## Pre-compute upstream quantities (single vectorised chemistry call)
    eval_gamma_star = chemistry["eval_gamma_star"]
    gamma_1 = eval_gamma_star(rho_1_vec, e_1_vec)

    e1_over_w1_sq = e_1_vec ./ (w_1_vec .^ 2)
    const_1 = (gamma_1 .- 1.0) .* e1_over_w1_sq .+ 1.0
    const_2 = gamma_1 .* e1_over_w1_sq .+ 0.5

    ## Build initial guess
    x0 = create_initial_guess(N, gamma_1, w_1_vec, e_1_vec)

    ## Define residual function for NLsolve
    function equations_vectorized!(F, x)
        # Reshape x to [N x 2] matrix: [p1 q1; p2 q2; ... ; pN qN]
        vars = reshape(x, 2, N)
        p_vals = vars[1, :]
        q_vals = vars[2, :]

        # Post-shock state for all systems
        rho_2_vals = p_vals .* rho_1_vec
        e_2_vals = q_vals .* (w_1_vec .^ 2)

        # Single vectorised chemistry call -- major speed-up
        gamma_2 = eval_gamma_star(rho_2_vals, e_2_vals)

        # Residuals for conservation equations
        eq1 = const_1 .- ((gamma_2 .- 1.0) .* p_vals .* q_vals .+ 1.0 ./ p_vals)
        eq2 = const_2 .- (gamma_2 .* q_vals .+ 1.0 ./ (2.0 .* p_vals .^ 2))

        # Interleave: [eq1_1, eq2_1, eq1_2, eq2_2, ..., eq1_N, eq2_N]
        for i in 1:N
            F[2*(i-1)+1] = eq1[i]
            F[2*(i-1)+2] = eq2[i]
        end
    end

    ## Configure and run NLsolve
    result = nlsolve(equations_vectorized!, x0,
                     method=:trust_region,
                     iterations=100,
                     ftol=1e-8,
                     xtol=1e-8,
                     show_trace=false)

    x_sol = result.zero

    ## Extract and compute post-shock state
    vars_sol = reshape(x_sol, 2, N)
    p_sol = vars_sol[1, :]
    q_sol = vars_sol[2, :]

    rho_2 = p_sol .* rho_1_vec
    e_2 = q_sol .* (w_1_vec .^ 2)
    w_2 = w_1_vec ./ p_sol

    return rho_2, e_2, w_2
end

# ========================================================================
#  Helper: physics-based initial guess
# ========================================================================
"""
    create_initial_guess(N, gamma_1, w_1, e_1)

Compute a physics-based starting point for NLsolve using ideal-gas
normal-shock relations.
"""
function create_initial_guess(N::Int, gamma_1::Vector{Float64},
                              w_1::Vector{Float64}, e_1::Vector{Float64})
    # Approximate Mach number squared
    M1_sq_approx = w_1 .^ 2 ./ (gamma_1 .* (gamma_1 .- 1.0) .* e_1)

    # Density-ratio guess from normal-shock relation
    p_guess = (gamma_1 .+ 1.0) .* M1_sq_approx ./ ((gamma_1 .- 1.0) .* M1_sq_approx .+ 2.0)

    # Energy-ratio guess
    q_guess = (0.5 .+ e_1 ./ w_1 .^ 2) ./ gamma_1

    if N == 1
        x0 = [p_guess[1]; q_guess[1]]
    else
        # Interleave: [p1, q1, p2, q2, ..., pN, qN]
        x0 = zeros(2 * N)
        for i in 1:N
            x0[2*(i-1)+1] = p_guess[i]
            x0[2*(i-1)+2] = q_guess[i]
        end
    end

    return x0
end
