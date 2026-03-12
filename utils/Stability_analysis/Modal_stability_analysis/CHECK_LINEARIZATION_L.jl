using LinearAlgebra
using Printf

"""
    CHECK_LINEARIZATION_L(A, s, chemistry)

Verify the linearized Jacobian A against finite differences.

Validates the Jacobian matrix A by comparing A*q (the linearized response)
against a finite-difference approximation obtained from the nonlinear operator:
    dq = [f(x + Dx) - f(x)] / Dx
Random perturbations are applied to all conservative variables
(rho, rho_u, rho_v, rho_E) and optionally to the shock position.

# Arguments
- `A`: Linearized Jacobian matrix (sparse).
- `s`: Solution Dict{String,Any} containing the base flow, grid, and solver parameters.
- `chemistry`: Chemistry model Dict{String,Any} for thermodynamic evaluations.

# Returns
- Nothing. Prints the relative linearization error to the console.

# Notes
- The perturbation magnitude is controlled by `s["stability_analysis"]["perturbation_magnitude"]`.
- When `s["stability_analysis"]["perturb_shock"]` is true, the shock
  position is also perturbed in the normal direction.
- The nonlinear operator is evaluated via NON_LINEAR_DYNAMICS_NO_DISCONTINUITY.
- A small value of the reported error (e.g., < 1e-3) indicates
  that the linearization is consistent with the nonlinear operator.

Part of: Hypersonics Stability Julia Solver - Stability Analysis / Linear Stability Module
"""
function CHECK_LINEARIZATION_L(A, s::Dict{String,Any}, chemistry::Dict{String,Any})
    ## Initialize perturbation and response vectors
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    if s["stability_analysis"]["perturb_shock"]
        q = zeros(4 * Nchi * Neta + Nchi, 1)
        dq = zeros(4 * Nchi * Neta + Nchi, 1)
    else
        q = zeros(4 * Nchi * Neta, 1)
        dq = zeros(4 * Nchi * Neta, 1)
    end

    solution_perturbed = deepcopy(s)

    ## Perturb rho
    matrix_perturbation = randn(Nchi, Neta)
    if s["shock"]["enabled"]
        solution_perturbed["var"]["rho"][2:end-1, 2:end-1] .= s["var"]["rho"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[1:Nchi*Neta, 1] .= reshape(s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    else
        solution_perturbed["var"]["rho"][2:end-1, 2:end-1] .= s["var"]["rho"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[1:Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    end

    ## Perturb rho_u
    matrix_perturbation = randn(Nchi, Neta)
    if s["shock"]["enabled"]
        solution_perturbed["var"]["rho_u"][2:end-1, 2:end-1] .= s["var"]["rho_u"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[Nchi*Neta+1:2*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    else
        solution_perturbed["var"]["rho_u"][2:end-1, 2:end-1] .= s["var"]["rho_u"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[Nchi*Neta+1:2*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    end

    ## Perturb rho_v
    matrix_perturbation = randn(Nchi, Neta)
    if s["shock"]["enabled"]
        solution_perturbed["var"]["rho_v"][2:end-1, 2:end-1] .= s["var"]["rho_v"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[2*Nchi*Neta+1:3*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    else
        solution_perturbed["var"]["rho_v"][2:end-1, 2:end-1] .= s["var"]["rho_v"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[2*Nchi*Neta+1:3*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    end

    ## Perturb rho_E
    matrix_perturbation = randn(Nchi, Neta)
    if s["shock"]["enabled"] && s["shock"]["feedback"]
        solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    elseif s["shock"]["enabled"] && !s["shock"]["feedback"]
        solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ s["shock"]["flow_cells"] .* matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    else
        solution_perturbed["var"]["rho_E"][2:end-1, 2:end-1] .= s["var"]["rho_E"][2:end-1, 2:end-1] .+ matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"]
        q[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape(matrix_perturbation .* s["stability_analysis"]["perturbation_magnitude"], Nchi*Neta)
    end

    ## Perturb shock position (if enabled)
    if s["stability_analysis"]["perturb_shock"]
        vector_perturbation = randn(Nchi, 1)
        ang = s["shock"]["beta"] .- pi/2 .- atan.(s["freestream"]["rho_v_0"], s["freestream"]["rho_u_0"])

        solution_perturbed["shock"]["points_x"] .= solution_perturbed["shock"]["points_x"] .+ s["stability_analysis"]["perturbation_magnitude"] .* cos.(ang) .* vector_perturbation
        solution_perturbed["shock"]["points_y"] .= solution_perturbed["shock"]["points_y"] .+ s["stability_analysis"]["perturbation_magnitude"] .* sin.(ang) .* vector_perturbation
        q[4*Nchi*Neta+1:4*Nchi*Neta+Nchi, 1] .= s["stability_analysis"]["perturbation_magnitude"] .* vector_perturbation
    end

    ## Evaluate nonlinear operator on base and perturbed solutions
    s_eval = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry)
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)

    ## Construct finite-difference response vector dq
    if s_eval["shock"]["enabled"]
        dq[0*Nchi*Neta+1:1*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho"] .- s_eval["flux"]["rho"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
        dq[1*Nchi*Neta+1:2*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_u"] .- s_eval["flux"]["rho_u"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
        dq[2*Nchi*Neta+1:3*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_v"] .- s_eval["flux"]["rho_v"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
        dq[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_E"] .- s_eval["flux"]["rho_E"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
    else
        dq[0*Nchi*Neta+1:1*Nchi*Neta, 1] .= reshape(solution_perturbed["flux"]["rho"] .- s_eval["flux"]["rho"], Nchi*Neta)
        dq[1*Nchi*Neta+1:2*Nchi*Neta, 1] .= reshape(solution_perturbed["flux"]["rho_u"] .- s_eval["flux"]["rho_u"], Nchi*Neta)
        dq[2*Nchi*Neta+1:3*Nchi*Neta, 1] .= reshape(solution_perturbed["flux"]["rho_v"] .- s_eval["flux"]["rho_v"], Nchi*Neta)
        dq[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape(solution_perturbed["flux"]["rho_E"] .- s_eval["flux"]["rho_E"], Nchi*Neta)
    end

    if s_eval["stability_analysis"]["perturb_shock"]
        dq[4*Nchi*Neta+1:4*Nchi*Neta+Nchi, 1] .= solution_perturbed["shock"]["relative_increase_velocity"] .- s_eval["shock"]["relative_increase_velocity"]
    end

    ## Compare linearized response A*q against finite-difference dq
    dq_lin = A * q

    temp = (dq_lin .- dq) ./ (abs.(dq_lin) .+ 1e-8)
    error_val = norm((dq_lin .- dq) ./ (abs.(dq_lin) .+ 1e-8)) / sqrt(9 * 4 * 4 * Nchi * Neta)
    @printf("Error of linearization: %e . \n", error_val)

    return nothing
end
