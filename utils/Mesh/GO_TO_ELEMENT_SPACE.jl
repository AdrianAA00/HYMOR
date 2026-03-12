# GO_TO_ELEMENT_SPACE  Convert physical (x,y) coordinates to element space [0,1]x[0,1].
#
#   (chi, eta) = GO_TO_ELEMENT_SPACE(x, y, s)
#
#   Transforms physical-space Cartesian coordinates (x, y) into the
#   element-space parametric coordinates (chi, eta) on the unit square
#   [0,1]x[0,1]. This is the inverse of GO_TO_PHYSICAL_SPACE.
#
#   Inputs:
#       x   - Array of physical x-coordinates (any size)
#       y   - Array of physical y-coordinates (same size as x)
#       s   - Dict containing curvilinear_mapping parameters
#
#   Outputs:
#       chi - Element-space streamwise coordinate in [0,1]
#       eta - Element-space wall-normal coordinate in [0,1]
#
#   Notes:
#       - For "circle", standard polar coordinates are converted to [0,1].
#       - For "MSL"/"blunt_cone", the piecewise geometry inverse is
#         applied, then scaled to [0,1].
#       - For "channel", "wedge", the inverse accounts for distance_max.
#       - For "lid_driven_cavity", Newton iteration inverts the nonlinear
#         LDC mapping.
#
# Part of: Hypersonics Stability Julia Solver - Mesh Module

function GO_TO_ELEMENT_SPACE(x::AbstractArray{Float64}, y::AbstractArray{Float64}, s::Dict{String,Any})

    boundary_type = s["curvilinear_mapping"]["boundary_type"]

    if boundary_type == "circle"
        ## Circle (polar coordinates -> [0,1])
        #  chi=0 at left boundary (angle_max), chi=1 at right boundary (angle_min)
        angle_min = -s["curvilinear_mapping"]["circle_angle_extra"]
        angle_max = pi/2

        angle = mod.(atan.(-x, y) .+ pi/2, 2*pi) .- pi/2
        chi = (angle_max .- angle) ./ (angle_max - angle_min)

        r = sqrt.(x.^2 .+ y.^2)
        distance_max = s["curvilinear_mapping"]["dRs"] .+
            (s["curvilinear_mapping"]["dRe"] - s["curvilinear_mapping"]["dRs"]) .* chi.^2
        eta = (r .- s["curvilinear_mapping"]["R"]) ./ distance_max

    elseif boundary_type == "MSL" || boundary_type == "blunt_cone"
        ## MSL / blunt cone (piecewise geometry -> [0,1])
        #  chi=0 at left boundary (chi_arc=0.5), chi=1 at right boundary (chi_arc=0)
        size_chi = size(x, 1)
        eta_phys = x .* 0
        chi_arc = x .* 0
        s_Straight = ((s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"]) +
            s["curvilinear_mapping"]["R"] * cos(s["curvilinear_mapping"]["theta"])) / sin(s["curvilinear_mapping"]["theta"])
        s_Curve = 2 * s["curvilinear_mapping"]["theta"] * s["curvilinear_mapping"]["R"]
        s_tot = s_Curve + 2 * s_Straight

        # Physical to body-aligned frame: x_body = y, y_body = -x
        x_body = y
        y_body = -x

        for i = 1:size_chi
            tang = (y_body[i,1] - (s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"])) / x_body[i,1]

            if (tang < tan(pi/2 - s["curvilinear_mapping"]["theta"])) && (x_body[i,1] > 0)
                # Right straight segment
                eta_phys[i,:] = y_body[i,:] .* cos(s["curvilinear_mapping"]["theta"]) .+
                    x_body[i,:] .* sin(s["curvilinear_mapping"]["theta"]) .-
                    s["curvilinear_mapping"]["R"] .-
                    (s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"]) * cos(s["curvilinear_mapping"]["theta"])
                lambda = (x_body[i,:] ./ cos(s["curvilinear_mapping"]["theta"]) .-
                    (eta_phys[i,:] .* sin(s["curvilinear_mapping"]["theta"]) .+
                    s["curvilinear_mapping"]["R"] * sin(s["curvilinear_mapping"]["theta"])) ./
                    cos(s["curvilinear_mapping"]["theta"])) ./ s_tot
                chi_arc[i,:] = s_Straight / s_tot .- lambda
            elseif (-tan(pi/2 - s["curvilinear_mapping"]["theta"]) < tang) && (x_body[i,1] < 0)
                # Left straight segment
                eta_phys[i,:] = y_body[i,:] .* cos(s["curvilinear_mapping"]["theta"]) .-
                    x_body[i,:] .* sin(s["curvilinear_mapping"]["theta"]) .-
                    s["curvilinear_mapping"]["R"] .-
                    (s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"]) * cos(s["curvilinear_mapping"]["theta"])
                lambda = (-x_body[i,:] ./ cos(s["curvilinear_mapping"]["theta"]) .-
                    (eta_phys[i,:] .* sin(s["curvilinear_mapping"]["theta"]) .+
                    s["curvilinear_mapping"]["R"] * sin(s["curvilinear_mapping"]["theta"])) ./
                    cos(s["curvilinear_mapping"]["theta"])) ./ s_tot
                chi_arc[i,:] = (s_Straight + s_Curve) / s_tot .+ lambda
            else
                # Circular nose cap
                angle = atan.((y_body[i,:] .+ s["curvilinear_mapping"]["R"] .- s["curvilinear_mapping"]["L"]), x_body[i,:])
                eta_phys[i,:] = x_body[i,:] ./ cos.(angle) .- s["curvilinear_mapping"]["R"]
                chi_arc[i,:] = (angle .- pi/2 .+ s["curvilinear_mapping"]["theta"]) ./
                    (2 * s["curvilinear_mapping"]["theta"] * s_tot / s_Curve) .+ s_Straight / s_tot
            end
        end

        # Convert from arc-length coordinates to [0,1] (inverted: chi_arc=0.5 -> chi=0, chi_arc=0 -> chi=1)
        chi = 1 .- 2 .* chi_arc

        # Compute distance_max and normalize eta
        if s["curvilinear_mapping"]["boundary_type"] == "blunt_cone"
            s_Str_sh = ((s["curvilinear_mapping"]["L"] - s["curvilinear_mapping"]["R"]) +
                s["curvilinear_mapping"]["R"] * cos(s["shock"]["initial_beta"])) / sin(s["shock"]["initial_beta"])
            s_Crv_sh = s["shock"]["initial_beta"] * s["curvilinear_mapping"]["R"]
            s_tot_sh = 2 * s_Crv_sh + 2 * s_Str_sh
            slope = sin(s["curvilinear_mapping"]["theta"] - s["shock"]["initial_beta"])

            distance_max = zeros(size(chi))
            for i = 1:size_chi
                if chi_arc[i,1] < s_Str_sh / s_tot_sh
                    distance_max[i,:] .= s["curvilinear_mapping"]["dRs"] +
                        (s_Str_sh - chi_arc[i,1] * s_tot_sh) * slope
                else
                    distance_max[i,:] .= s["curvilinear_mapping"]["dRs"]
                end
            end
        else
            distance_max = s["curvilinear_mapping"]["dRs"] .+
                (s["curvilinear_mapping"]["dRe"] - s["curvilinear_mapping"]["dRs"]) .* chi.^2
        end

        eta = eta_phys ./ distance_max

    elseif boundary_type == "channel"
        ## Channel (inverse with channel angle)
        chi = (x .* cos(s["curvilinear_mapping"]["channel_angle"]) .+ y .* sin(s["curvilinear_mapping"]["channel_angle"])) ./ s["curvilinear_mapping"]["Lx"]
        eta_r = (-x .* sin(s["curvilinear_mapping"]["channel_angle"]) .+ y .* cos(s["curvilinear_mapping"]["channel_angle"])) ./ s["curvilinear_mapping"]["Ly"]

        # Invert wall-normal refinement (p = eta_refinement_power)
        p = s["curvilinear_mapping"]["eta_refinement_power"]
        eta = zeros(size(eta_r))
        mask = eta_r .<= 0.5
        eta[mask]  = (eta_r[mask]  ./ 2^(p-1)).^(1/p)
        eta[.!mask] = 1 .- ((1 .- eta_r[.!mask]) ./ 2^(p-1)).^(1/p)

    elseif boundary_type == "lid_driven_cavity"
        ## Lid-driven cavity (Newton inversion of nonlinear LDC mapping)
        chi_phys = 1 .- y
        eta_phys = -x

        r = s["curvilinear_mapping"]["refinement_ratio_LDC"]
        f_max = 1 + 2*(r-1)/3

        # Invert F(s)/f_max = target via Newton iteration
        chi = ldc_inverse(chi_phys, r, f_max)
        eta = ldc_inverse(eta_phys, r, f_max)

    elseif boundary_type == "wedge"
        ## Wedge / Linear Interaction Analysis
        chi = 1 .- y ./ s["curvilinear_mapping"]["Lx"]
        distance_max = s["curvilinear_mapping"]["dRs"] .+
            (s["curvilinear_mapping"]["dRe"] - s["curvilinear_mapping"]["dRs"]) .* (1 .- chi)
        eta = -x ./ distance_max

    else
        error("GO_TO_ELEMENT_SPACE type not defined")
    end

    return chi, eta
end


# ldc_inverse  Invert the LDC nonlinear mapping via Newton iteration.
#
#   Given target = F(s)/f_max, solve for s where
#   F(s) = s + 4*(r-1)*(s^2/2 - s^3/3).

function ldc_inverse(target::AbstractArray{Float64}, r::Float64, f_max::Float64)

    target_scaled = target .* f_max
    s_inv = copy(target)  # Initial guess (exact when r=1)
    F  = similar(s_inv)
    Fp = similar(s_inv)

    coeff = 4*(r-1)
    for iter = 1:10
        @. F  = s_inv + coeff * (s_inv^2 / 2 - s_inv^3 / 3)
        @. Fp = 1 + coeff * s_inv * (1 - s_inv)
        @. s_inv = s_inv - (F - target_scaled) / Fp
    end

    return s_inv
end
