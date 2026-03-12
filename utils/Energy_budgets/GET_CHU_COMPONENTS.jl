function GET_CHU_COMPONENTS(
        rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
        pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
    # GET_CHU_COMPONENTS  Transform conservative perturbations to Chu energy variables.
    #
    #   (pert_p, pert_u, pert_v, pert_entropy) = GET_CHU_COMPONENTS(
    #       rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
    #       pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
    #
    #   Converts perturbation quantities from the conservative variable set
    #   (rho', rho_u', rho_v', rho_E') to the primitive / Chu energy variable
    #   set (p', u', v', s'/R) used in the Chu energy norm decomposition.
    #
    #   Inputs:
    #       rho_0        - Base-flow density.
    #       rho_u_0      - Base-flow x-momentum (rho * u).
    #       rho_v_0      - Base-flow y-momentum (rho * v).
    #       p_0          - Base-flow pressure.
    #       a_0          - Base-flow speed of sound (unused in current formulation
    #                      but retained for interface consistency with GET_ENERGY_NORM).
    #       gamma_star_0 - Base-flow effective ratio of specific heats (gamma*).
    #       pert_rho     - Perturbation density (rho').
    #       pert_rho_u   - Perturbation x-momentum (rho_u').
    #       pert_rho_v   - Perturbation y-momentum (rho_v').
    #       pert_rho_E   - Perturbation total energy (rho_E').
    #
    #   Outputs:
    #       pert_p       - Perturbation pressure.
    #       pert_u       - Perturbation x-velocity.
    #       pert_v       - Perturbation y-velocity.
    #       pert_entropy - Non-dimensional entropy perturbation (s'/R).
    #
    #   Notes:
    #       - All input arrays must be the same size; the function operates
    #         element-wise and supports arbitrary array dimensions.
    #       - The entropy output is normalised by the gas constant R so that
    #         it is a dimensionless quantity (s'/R). To recover s', multiply
    #         by R_0 = (gamma*-1) * cv*.
    #       - The parameter a_0 (speed of sound) is accepted for interface
    #         consistency with GET_ENERGY_NORM but is not used in the
    #         transformation formulae.

    # Base-flow primitive velocities
    u_0 = rho_u_0 ./ rho_0
    v_0 = rho_v_0 ./ rho_0

    # Velocity perturbations
    pert_u = (pert_rho_u .- u_0 .* pert_rho) ./ rho_0
    pert_v = (pert_rho_v .- v_0 .* pert_rho) ./ rho_0

    # Pressure perturbation from the energy equation
    t1 = -u_0 .* pert_rho_u
    t2 = -v_0 .* pert_rho_v
    t3 = (u_0.^2 .+ v_0.^2) .* pert_rho ./ 2
    pert_p = (gamma_star_0 .- 1) .* (pert_rho_E .+ t1 .+ t2 .+ t3)

    # Non-dimensional entropy perturbation (s'/R)
    pert_entropy = pert_p ./ p_0 ./ (gamma_star_0 .- 1) .-
                   pert_rho ./ rho_0 .* gamma_star_0 ./ (gamma_star_0 .- 1)

    return (pert_p, pert_u, pert_v, pert_entropy)
end
