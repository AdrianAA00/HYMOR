function GET_ENERGY_NORM(rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
                         pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
    # GET_ENERGY_NORM  Compute the Chu energy norm components for a perturbation field.
    #
    #   E = GET_ENERGY_NORM(rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
    #                       pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)
    #
    #   Evaluates the three components (acoustic, kinetic, entropic) of the Chu
    #   perturbation energy norm, which provides a positive-definite measure of
    #   disturbance energy in compressible flows.
    #
    #   Inputs:
    #       rho_0        - Base-flow density.
    #       rho_u_0      - Base-flow x-momentum (rho * u).
    #       rho_v_0      - Base-flow y-momentum (rho * v).
    #       p_0          - Base-flow pressure.
    #       a_0          - Base-flow speed of sound.
    #       gamma_star_0 - Base-flow effective ratio of specific heats (gamma*).
    #       pert_rho     - Perturbation density (rho').
    #       pert_rho_u   - Perturbation x-momentum (rho_u').
    #       pert_rho_v   - Perturbation y-momentum (rho_v').
    #       pert_rho_E   - Perturbation total energy (rho_E').
    #
    #   Outputs:
    #       E - (Dict{String,Any}) Chu energy norm components (same size as inputs):
    #           "acoustic" - Acoustic energy
    #           "kinetic"  - Kinetic energy
    #           "entropic" - Entropic energy
    #
    #   Notes:
    #       - Conservative perturbations are first transformed to Chu variables
    #         (p', u', v', s'/R) via GET_CHU_COMPONENTS.
    #       - Only the real parts of the perturbations are used in the quadratic
    #         energy expressions.

    ## Transform to Chu variables
    (pert_p, pert_u, pert_v, pert_entropy) = GET_CHU_COMPONENTS(
        rho_0, rho_u_0, rho_v_0, p_0, a_0, gamma_star_0,
        pert_rho, pert_rho_u, pert_rho_v, pert_rho_E)

    ## Evaluate energy norm components
    E = Dict{String,Any}()

    # Acoustic energy: rho_0 * a_0^2 / (2 * (gamma* * p_0)^2) * p'^2
    E["acoustic"] = rho_0 .* a_0.^2 ./ (2 .* (gamma_star_0 .* p_0).^2) .* real.(pert_p).^2

    # Kinetic energy: rho_0 / 2 * (u'^2 + v'^2)
    E["kinetic"] = rho_0 ./ 2 .* real.(pert_u).^2 .+ rho_0 ./ 2 .* real.(pert_v).^2

    # Entropic energy: (gamma* - 1) * p_0 / (2 * gamma*) * (s'/R)^2
    E["entropic"] = (gamma_star_0 .- 1) .* p_0 ./ (2 .* gamma_star_0) .* real.(pert_entropy).^2

    return E
end
