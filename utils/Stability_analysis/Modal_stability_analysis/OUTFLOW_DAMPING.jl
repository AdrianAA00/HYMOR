# OUTFLOW_DAMPING  Damping (sponge) layer of the linear operator at outflow boundaries.
#
#   (D_damp, sigma) = OUTFLOW_DAMPING(s)
#
#   Builds the diagonal operator of a damping layer, d q'/dt = ... - sigma q',
#   in a thin layer of cells next to every boundary whose condition is of
#   outflow type ("outflow_supersonic", "outflow_subsonic", "outflow_NRCBC",
#   "outflow_supersonic_1st"). Disturbances reaching such a boundary are
#   damped before they can reflect, so the slowly decaying (or spuriously
#   growing) modes that live at a truncated outflow do not dominate the
#   spectrum. The base flow is not modified: the layer only acts on
#   perturbations, i.e. it is subtracted from the linearized operator L.
#   Cells outside the flow (beyond the fitted shock), which have no
#   dynamics, are damped at the same strength so that their lambda = 0
#   eigenspace does not crowd the leading spectrum.
#
#   The layer spans a fraction of the cells normal to the boundary (so its
#   width in chi/eta is the same on every grid) and sigma ramps smoothly
#   from 0 at its inner edge to the strength at the boundary:
#       sigma = strength * sin(pi/2 * xi)^2,   xi = 0 (inner edge) ... 1 (boundary)
#
#   Inputs:
#       s - Solution Dict with:
#           s["stability_analysis"]["outflow_damping"]["enabled"]  - (Bool)
#           s["stability_analysis"]["outflow_damping"]["fraction"] - layer width, fraction of the cells
#           s["stability_analysis"]["outflow_damping"]["strength"] - sigma at the boundary (1/time)
#           s["boundary_conditions"], s["mesh"]["Nchi"], s["mesh"]["Neta"]
#
#   Outputs:
#       D_damp - Sparse diagonal matrix of size of the linear state
#                (rho, rho_u, rho_v, rho_E blocks, then shock DoFs, which are
#                not damped). Zero matrix if the layer is disabled or no
#                boundary is of outflow type.
#       sigma  - (Nchi x Neta) damping rate of each cell (layer and cells
#                outside the flow).
#
# Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Linear Stability Module

function OUTFLOW_DAMPING(s::Dict{String,Any})
    Nchi = s["mesh"]["Nchi"]
    Neta = s["mesh"]["Neta"]
    sigma = zeros(Nchi, Neta)

    damping = s["stability_analysis"]["outflow_damping"]
    outflow_types = ("outflow_supersonic", "outflow_subsonic", "outflow_NRCBC", "outflow_supersonic_1st")
    if damping["enabled"]
        for b in ("chi0", "chi1", "eta0", "eta1")
            if !(s["boundary_conditions"]["boundary_" * b]["name"] in outflow_types)
                continue
            end
            n_normal = startswith(b, "chi") ? Nchi : Neta
            n_layer = max(2, ceil(Int, damping["fraction"] * n_normal))
            xi = ((1:n_layer) .- 0.5) ./ n_layer                # 0 at the inner edge, 1 at the boundary
            ramp = damping["strength"] .* sin.(pi / 2 .* xi) .^ 2
            if b == "chi0"
                sigma[1:n_layer, :] .= max.(sigma[1:n_layer, :], reverse(ramp) * ones(1, Neta))
            elseif b == "chi1"
                sigma[end-n_layer+1:end, :] .= max.(sigma[end-n_layer+1:end, :], ramp * ones(1, Neta))
            elseif b == "eta0"
                sigma[:, 1:n_layer] .= max.(sigma[:, 1:n_layer], ones(Nchi) * reverse(ramp)')
            else
                sigma[:, end-n_layer+1:end] .= max.(sigma[:, end-n_layer+1:end], ones(Nchi) * ramp')
            end
        end

        # Cells outside the flow (beyond the fitted shock) carry no dynamics:
        # their rows and columns of L vanish, so they add a large lambda = 0
        # eigenspace that crowds the leading spectrum once the outflow modes
        # are damped. Damping them as well moves that eigenspace to -strength.
        flow_cells = s["shock"]["flow_cells"]
        if s["shock"]["enabled"] && flow_cells isa AbstractMatrix && size(flow_cells) == (Nchi, Neta)
            outside = flow_cells .== 0
            sigma[outside] .= max.(sigma[outside], damping["strength"])
        end
    end

    n_shock = s["stability_analysis"]["perturb_shock"] ? Nchi : 0
    diagonal = vcat(repeat(vec(sigma), 4), zeros(n_shock))
    D_damp = spdiagm(0 => diagonal)
    return (D_damp, sigma)
end
