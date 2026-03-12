using LinearAlgebra
using Printf

"""
    CHECK_LINEARIZATION_B(B, s, chemistry)

Verify the linearized upstream boundary operator B by finite differences.

Validates the Jacobian matrix B (upstream boundary linearization) by
comparing B*q against a finite-difference approximation computed from
the full nonlinear dynamics. Prints the relative error norm.

# Arguments
- `B`: Linearized upstream boundary operator (sparse matrix).
- `s`: Solution Dict{String,Any} with base flow, mesh, and parameters.
- `chemistry`: Chemistry model Dict{String,Any}.

# Returns
- Nothing. Prints the linearization error to the console.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function CHECK_LINEARIZATION_B(B, s::Dict{String,Any}, chemistry::Dict{String,Any})
    ## Work on a local copy to match MATLAB pass-by-value semantics and
    ## prevent rho_0_upstream_p / freestream.*_p from leaking to the caller.
    s = deepcopy(s)

    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]

    ## Set uniform upstream base state
    s["freestream"]["rho_0_p"] = s["freestream"]["rho_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_u_0_p"] = s["freestream"]["rho_u_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_v_0_p"] = s["freestream"]["rho_v_0"] .* ones(Nchi, 1)
    s["freestream"]["rho_E_0_p"] = s["freestream"]["rho_E_0"] .* ones(Nchi, 1)
    s["rho_0_upstream_p"] = true

    q = zeros(4 * Nchi, 1)
    N_dof = s["stability_analysis"]["perturb_shock"] ? 4 * Nchi * Neta + Nchi : 4 * Nchi * Neta
    dq = zeros(N_dof, 1)

    ## Compute nonlinear finite-difference perturbation
    solution_perturbed = deepcopy(s)

    # rho perturbation
    p_pert = randn(Nchi, 1)
    solution_perturbed["freestream"]["rho_0_p"] = s["freestream"]["rho_0_p"] .+ p_pert .* s["stability_analysis"]["perturbation_magnitude"]
    q[1:Nchi, 1] .= p_pert .* s["stability_analysis"]["perturbation_magnitude"]

    # rho_u perturbation
    p_pert = randn(Nchi, 1)
    solution_perturbed["freestream"]["rho_u_0_p"] = s["freestream"]["rho_u_0_p"] .+ p_pert .* s["stability_analysis"]["perturbation_magnitude"]
    q[1+Nchi:2*Nchi, 1] .= p_pert .* s["stability_analysis"]["perturbation_magnitude"]

    # rho_v perturbation
    p_pert = randn(Nchi, 1)
    solution_perturbed["freestream"]["rho_v_0_p"] = s["freestream"]["rho_v_0_p"] .+ p_pert .* s["stability_analysis"]["perturbation_magnitude"]
    q[1+2*Nchi:3*Nchi, 1] .= p_pert .* s["stability_analysis"]["perturbation_magnitude"]

    # rho_E perturbation
    p_pert = randn(Nchi, 1)
    solution_perturbed["freestream"]["rho_E_0_p"] = s["freestream"]["rho_E_0_p"] .+ p_pert .* s["stability_analysis"]["perturbation_magnitude"]
    q[1+3*Nchi:4*Nchi, 1] .= p_pert .* s["stability_analysis"]["perturbation_magnitude"]

    ## Evaluate nonlinear dynamics for base and perturbed states
    s_eval = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(s, chemistry)
    solution_perturbed = NON_LINEAR_DYNAMICS_NO_DISCONTINUITY(solution_perturbed, chemistry)

    ## Build finite-difference Jacobian-vector product
    dq[0*Nchi*Neta+1:1*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho"] .- s_eval["flux"]["rho"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
    dq[1*Nchi*Neta+1:2*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_u"] .- s_eval["flux"]["rho_u"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
    dq[2*Nchi*Neta+1:3*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_v"] .- s_eval["flux"]["rho_v"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)
    dq[3*Nchi*Neta+1:4*Nchi*Neta, 1] .= reshape((solution_perturbed["flux"]["rho_E"] .- s_eval["flux"]["rho_E"]) .* s_eval["shock"]["flow_cells"], Nchi*Neta)

    if s_eval["stability_analysis"]["perturb_shock"]
        dq[4*Nchi*Neta+1:4*Nchi*Neta+Nchi, 1] .= solution_perturbed["shock"]["relative_increase_velocity"] .- s_eval["shock"]["relative_increase_velocity"]
    end

    ## Compare linearized vs finite-difference result
    dq_lin = B * q

    ## Compute and report error
    temp = (dq_lin .- dq) ./ (abs.(dq_lin) .+ 1e-8)
    (a_val, b_val) = findmax(temp)
    error_val = norm((dq_lin .- dq) ./ (abs.(dq_lin) .+ 1e-8)) / sqrt(3 * 4 * 4 * Nchi)
    @printf("Error of linearization: %e . \n", error_val)

    return nothing
end
