function ENERGY_INFLOW(V, s, w_infty)
    # ENERGY_INFLOW  Compute inflow energy flux and its per-unit-volume density.
    #
    #   (dE_dt_inflow, dE_dt_V_inflow) = ENERGY_INFLOW(V, s, w_infty)
    #
    #   Evaluates the total perturbation energy flux entering the domain through
    #   the freestream boundary. The flux operator D is constructed from the
    #   freestream coupling matrices (R_ and M_infty_) and applied to the global
    #   perturbation vector V. Individual acoustic, kinetic, and entropic
    #   contributions are computed and reported. The result is also normalised
    #   by the total post-shock volume to yield an energy density metric.
    #
    #   Inputs:
    #       V        - (Vector) Global perturbation state vector.
    #       s        - (Dict) Solution structure with grid information and
    #                  freestream parameters.
    #       w_infty  - (Vector) Freestream modal frequency/weight vector.
    #
    #   Outputs:
    #       dE_dt_inflow   - (Float64) Total perturbation energy flux
    #                        through the inflow boundary (V' * D * V).
    #       dE_dt_V_inflow - (Float64) Energy flux normalised by the total
    #                        post-shock volume (energy density rate).

    ## Construct inflow flux operators
    R_ = CONSTRUCT_R_(s, w_infty)
    T = 1  # Flux-of-energy mode
    scaling_non_temporal = false

    # Total energy flux operator
    norms_total = [1, 1, 1, 1]
    M_infty_total = CONSTRUCT_M_INFTY_(s, norms_total, T, w_infty, scaling_non_temporal)
    D_ = R_' * M_infty_total * R_

    # Acoustic-only operator
    norms_p = [1, 0, 0, 0]
    M_infty_p = CONSTRUCT_M_INFTY_(s, norms_p, T, w_infty, scaling_non_temporal)
    D_p = R_' * M_infty_p * R_

    # Kinetic-only operator
    norms_u = [0, 1, 1, 0]
    M_infty_u = CONSTRUCT_M_INFTY_(s, norms_u, T, w_infty, scaling_non_temporal)
    D_u = R_' * M_infty_u * R_

    # Entropic-only operator
    norms_s = [0, 0, 0, 1]
    M_infty_s = CONSTRUCT_M_INFTY_(s, norms_s, T, w_infty, scaling_non_temporal)
    D_s = R_' * M_infty_s * R_

    ## Evaluate inflow energy fluxes
    dE_dt_inflow   = V' * D_  * V
    dE_dt_inflow_p = V' * D_p * V
    dE_dt_inflow_s = V' * D_s * V
    dE_dt_inflow_u = V' * D_u * V

    println("Freestream acoustic: " * string(abs(dE_dt_inflow_p) / abs(dE_dt_inflow)))
    println("Freestream entropic: " * string(abs(dE_dt_inflow_s) / abs(dE_dt_inflow)))
    println("Freestream kinetic: "  * string(abs(dE_dt_inflow_u) / abs(dE_dt_inflow)))

    ## Compute total post-shock volume (vectorized)
    volume     = s["mesh"]["volume"]
    flow_cells = s["shock"]["flow_cells"]
    total_post_shock_volume = sum(volume .* flow_cells)

    # Energy density: inflow flux per unit volume
    dE_dt_V_inflow = dE_dt_inflow / total_post_shock_volume

    return (dE_dt_inflow, dE_dt_V_inflow)
end
