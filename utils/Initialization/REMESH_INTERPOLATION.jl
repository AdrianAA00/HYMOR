# REMESH_INTERPOLATION  Transfer the flow fields of solution_old onto the
#   mesh of s by cubic interpolation in the (i, j) index space of the old
#   structured grid.
#
#   s = REMESH_INTERPOLATION(s, solution_old)
#
#   Inputs:
#       s            - (Dict) Solution whose extended mesh (x_Ext, y_Ext) is the target
#       solution_old - (Dict) Solution holding the fields on the old mesh
#
#   Outputs:
#       s - Same Dict with s["var"] rho, rho_u, rho_v, rho_E (and gamma_star,
#           cv_star out of chemical equilibrium) filled on the extended grid.
#
#   Notes:
#       - Each new cell centre is located in the old grid by inverting the
#         map (i, j) -> (x, y) of the old extended cell centres with Newton
#         iterations (started from the nearest old node). The fields are
#         then interpolated with tensor-product cubic convolution (Keys,
#         a = -1/2) on the uniform (i, j) lattice: local 4 x 4 stencil, C1,
#         third-order accurate.
#       - Interpolating along the grid lines keeps the interpolation error
#         smooth on stretched, curved, high-aspect-ratio grids. Delaunay-
#         linear interpolation in (x, y) instead picks triangle diagonals
#         that vary from cell to cell, and its error is amplified by the
#         discrete derivatives (larger residual jump after a remesh, mostly
#         in the boundary layer and at the stagnation region).
#       - With shock fitting, the old fields from the shocked cell up are
#         replaced by a linear continuation of the last two flow cells, so
#         no stencil straddles the shock or uses the shocked-cell closure.
#         Corner ghost cells (not set by the boundary conditions) are filled
#         by bilinear extrapolation.
#       - New points outside the old grid (ghost and freestream cells) are
#         extrapolated; they are overwritten afterwards by the boundary
#         conditions and the shock treatment.
#
# Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function REMESH_INTERPOLATION(s::Dict{String,Any}, solution_old::Dict{String,Any})

    ## Fields to transfer
    fields = ["rho", "rho_u", "rho_v", "rho_E"]
    if !s["chemistry"]["chemical_equilibrium"]
        append!(fields, ["gamma_star", "cv_star"])
    end

    ## Locate the new cell centres in the old index space
    xo = solution_old["mesh"]["x_Ext"]
    yo = solution_old["mesh"]["y_Ext"]
    (xi, eta) = REMESH_LOCATE_IN_OLD_GRID(xo, yo, s["mesh"]["x_Ext"], s["mesh"]["y_Ext"])

    ## Interpolate every field
    for f in fields
        F = REMESH_PREPARE_OLD_FIELD(solution_old["var"][f], solution_old)
        s["var"][f] = REMESH_KEYS_EVAL(REMESH_KEYS_PAD(F), xi, eta, false)[1]
    end

    return s
end


# Invert the cubic-convolution map (i, j) -> (x, y) of the old grid at each
# new point with damped Newton iterations, started from the nearest old
# node: nearest old grid line (matched on the first interior row), then
# nearest node along that line.
function REMESH_LOCATE_IN_OLD_GRID(xo::AbstractMatrix, yo::AbstractMatrix,
                                   xn::AbstractMatrix, yn::AbstractMatrix)
    (Ni, Nj) = size(xo)
    (Mi, Mj) = size(xn)

    xi  = zeros(Mi, Mj)
    eta = zeros(Mi, Mj)
    @inbounds for i in 1:Mi
        k = argmin([(xn[i, 2] - xo[m, 2])^2 + (yn[i, 2] - yo[m, 2])^2 for m in 1:Ni])
        for j in 1:Mj
            xi[i, j]  = k
            eta[i, j] = argmin([(xn[i, j] - xo[k, m])^2 + (yn[i, j] - yo[k, m])^2 for m in 1:Nj])
        end
    end

    Px = REMESH_KEYS_PAD(xo)
    Py = REMESH_KEYS_PAD(yo)
    last_step = fill(Inf, Mi, Mj)
    for it in 1:50
        (X, X_i, X_j) = REMESH_KEYS_EVAL(Px, xi, eta, true)
        (Y, Y_i, Y_j) = REMESH_KEYS_EVAL(Py, xi, eta, true)
        max_step = 0.0
        @inbounds for q in eachindex(xi)
            rx = xn[q] - X[q]
            ry = yn[q] - Y[q]
            det = X_i[q] * Y_j[q] - X_j[q] * Y_i[q]
            dxi  = ( Y_j[q] * rx - X_j[q] * ry) / det
            deta = (-Y_i[q] * rx + X_i[q] * ry) / det
            step = max(1.0, max(abs(dxi), abs(deta)))   # at most one cell per iteration
            xi[q]  += dxi / step
            eta[q] += deta / step
            last_step[q] = abs(dxi) + abs(deta)
            max_step = max(max_step, last_step[q])
        end
        if max_step < 1e-10
            break
        end
    end
    located = isfinite.(xi) .& isfinite.(eta) .& (last_step .<= 1e-6)
    if !all(located[2:end-1, 2:end-1])   # ghost cells are reset by the BCs
        error("REMESH_INTERPOLATION: new cell centres could not be located in the old grid.")
    end

    return xi, eta
end


# Replace the shocked cell and everything beyond it by a linear
# continuation of the last two flow cells, and fill the four corner ghost
# cells. The shocked cell holds the Rankine-Hugoniot reconstruction, not a
# smooth continuation of the flow field: interpolating through it gives a
# larger residual in the first row below the shock (recomputed anyway by
# UPDATE_SHOCK_BC on the new mesh).
function REMESH_PREPARE_OLD_FIELD(F_old::AbstractMatrix, solution_old::Dict{String,Any})
    F = Matrix{Float64}(copy(F_old))
    (Ni, Nj) = size(F)
    if solution_old["shock"]["enabled"]
        js = vec(solution_old["shock"]["cell_indices"][:, 1])
        c_col = vcat(js[1], js, js[end]) .+ 1   # shocked cell (extended indexing)
        @inbounds for i in 1:Ni
            m1 = Int(c_col[i]) - 1              # last flow cell below the shock
            for k in 1:(Nj - m1)
                F[i, m1 + k] = F[i, m1] + k * (F[i, m1] - F[i, m1 - 1])
            end
        end
    end
    F[1, 1]   = F[2, 1]      + F[1, 2]      - F[2, 2]
    F[Ni, 1]  = F[Ni-1, 1]   + F[Ni, 2]     - F[Ni-1, 2]
    F[1, Nj]  = F[2, Nj]     + F[1, Nj-1]   - F[2, Nj-1]
    F[Ni, Nj] = F[Ni-1, Nj]  + F[Ni, Nj-1]  - F[Ni-1, Nj-1]
    return F
end


# Pad the grid data by one node per side with Keys' boundary rule.
function REMESH_KEYS_PAD(F::AbstractMatrix)
    (Ni, Nj) = size(F)
    P = zeros(Ni + 2, Nj + 2)
    P[2:end-1, 2:end-1] .= F
    @views P[1, 2:end-1]   .= 3 .* F[1, :]   .- 3 .* F[2, :]     .+ F[3, :]
    @views P[end, 2:end-1] .= 3 .* F[end, :] .- 3 .* F[end-1, :] .+ F[end-2, :]
    @views P[:, 1]   .= 3 .* P[:, 2]     .- 3 .* P[:, 3]     .+ P[:, 4]
    @views P[:, end] .= 3 .* P[:, end-1] .- 3 .* P[:, end-2] .+ P[:, end-3]
    return P
end


# Cubic-convolution weights (Keys 1981, a = -1/2) and their derivatives for
# the stencil k-1..k+2 at fractional position t from node k.
@inline function REMESH_KEYS_KERNEL(t::Float64)
    t2 = t * t; t3 = t2 * t
    return ((-t3 + 2t2 - t) / 2, (3t3 - 5t2 + 2) / 2, (-3t3 + 4t2 + t) / 2, (t3 - t2) / 2)
end

@inline function REMESH_KEYS_KERNEL_DT(t::Float64)
    t2 = t * t
    return ((-3t2 + 4t - 1) / 2, (9t2 - 10t) / 2, (-9t2 + 8t + 1) / 2, (3t2 - 2t) / 2)
end


# Evaluate the interpolant of the padded data P (and its index derivatives)
# at fractional indices (xi, eta) of the unpadded grid. Queries outside the
# grid are extrapolated with the boundary cubic.
function REMESH_KEYS_EVAL(P::Matrix{Float64}, xi::AbstractArray, eta::AbstractArray, derivatives::Bool)
    Ni = size(P, 1) - 2
    Nj = size(P, 2) - 2
    v   = similar(xi, Float64)
    v_i = derivatives ? similar(xi, Float64) : v
    v_j = derivatives ? similar(xi, Float64) : v
    @inbounds for q in eachindex(xi)
        if !isfinite(xi[q]) || !isfinite(eta[q])
            v[q] = NaN
            if derivatives
                v_i[q] = NaN; v_j[q] = NaN
            end
            continue
        end
        ki = clamp(floor(Int, xi[q]), 1, Ni - 1);  ti = xi[q] - ki
        kj = clamp(floor(Int, eta[q]), 1, Nj - 1); tj = eta[q] - kj
        wi = REMESH_KEYS_KERNEL(ti)
        wj = REMESH_KEYS_KERNEL(tj)
        acc = 0.0
        for a in 1:4, b in 1:4
            acc += wi[a] * wj[b] * P[ki + a - 1, kj + b - 1]
        end
        v[q] = acc
        if derivatives
            di = REMESH_KEYS_KERNEL_DT(ti)
            dj = REMESH_KEYS_KERNEL_DT(tj)
            acc_i = 0.0; acc_j = 0.0
            for a in 1:4, b in 1:4
                Pab = P[ki + a - 1, kj + b - 1]
                acc_i += di[a] * wj[b] * Pab
                acc_j += wi[a] * dj[b] * Pab
            end
            v_i[q] = acc_i
            v_j[q] = acc_j
        end
    end
    return v, v_i, v_j
end
