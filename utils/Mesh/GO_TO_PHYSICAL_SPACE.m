function [x, y] = GO_TO_PHYSICAL_SPACE(chi, eta, s)
% GO_TO_PHYSICAL_SPACE  Convert element-space [0,1]x[0,1] coordinates to physical (x,y).
%
%   [x, y] = GO_TO_PHYSICAL_SPACE(chi, eta, s)
%
%   Transforms element-space parametric coordinates (chi, eta) on the
%   unit square [0,1]x[0,1] into physical-space Cartesian coordinates
%   (x, y) appropriate for the selected boundary type. chi parameterizes
%   the wall (streamwise direction) and eta parameterizes the wall-normal
%   direction (eta=0 at wall, eta=1 at outer/shock boundary).
%
%   Inputs:
%       chi - Streamwise parametric coordinate in [0,1]
%       eta - Wall-normal parametric coordinate in [0,1]
%       s   - Struct containing curvilinear_mapping parameters
%
%   Outputs:
%       x - Physical x-coordinates (same size as inputs)
%       y - Physical y-coordinates (same size as inputs)
%
%   Notes:
%       - For "circle", chi maps to angle and eta maps to radial distance.
%       - For "MSL"/"blunt_cone", chi maps to normalized arc length and
%         eta maps to wall-normal distance via a piecewise geometry.
%       - For "channel", "lid_driven_cavity", "wedge", chi maps to
%         streamwise length and eta maps to wall-normal distance.
%
% Part of: Hypersonics Stability MATLAB Solver - Mesh Module

    switch s.curvilinear_mapping.boundary_type

        %% Circle (chi -> angle, eta -> radial distance from R)
        %  chi=0 at left boundary (angle_max), chi=1 at right boundary (angle_min)
        case "circle"
            angle_min = -s.curvilinear_mapping.circle_angle_extra;
            angle_max = pi/2;
            angle = angle_max - chi * (angle_max - angle_min);

            distance_max = s.curvilinear_mapping.dRs ...
                + (s.curvilinear_mapping.dRe - s.curvilinear_mapping.dRs) * chi.^2;
            r = s.curvilinear_mapping.R + eta .* distance_max;

            x = -r .* sin(angle);
            y =  r .* cos(angle);

        %% MSL / blunt cone (piecewise geometry with wall-normal extrusion)
        %  chi=0 at left boundary (chi_arc=0.5), chi=1 at right boundary (chi_arc=0)
        case {"MSL", "blunt_cone"}
            chi_arc = (1 - chi) / 2;  % Map chi in [0,1] to arc-length fraction in [0.5, 0]

            % Compute wall-normal distance scaling
            if s.curvilinear_mapping.boundary_type == "blunt_cone"
                s_Str_sh = ((s.curvilinear_mapping.L - s.curvilinear_mapping.R) ...
                    + s.curvilinear_mapping.R * cos(s.shock.initial_beta)) / sin(s.shock.initial_beta);
                s_Crv_sh = s.shock.initial_beta * s.curvilinear_mapping.R;
                s_tot_sh = 2 * s_Crv_sh + 2 * s_Str_sh;
                slope = sin(s.curvilinear_mapping.theta - s.shock.initial_beta); % Extra slope when measured from shock initial position

                size_chi = size(chi, 1);
                distance_max = zeros(size(chi));
                for i = 1:size_chi
                    if chi_arc(i,1) < s_Str_sh / s_tot_sh
                        distance_max(i,:) = s.curvilinear_mapping.dRs ...
                            + (s_Str_sh - chi_arc(i,1) * s_tot_sh) * slope;
                    else
                        distance_max(i,:) = s.curvilinear_mapping.dRs;
                    end
                end
            else
                distance_max = s.curvilinear_mapping.dRs ...
                    + (s.curvilinear_mapping.dRe - s.curvilinear_mapping.dRs) * chi.^2;
            end

            eta_phys = eta .* distance_max;

            % Piecewise geometry parameters
            s_Straight = ((s.curvilinear_mapping.L - s.curvilinear_mapping.R) ...
                + s.curvilinear_mapping.R * cos(s.curvilinear_mapping.theta)) / sin(s.curvilinear_mapping.theta);
            s_Curve = 2 * s.curvilinear_mapping.theta * s.curvilinear_mapping.R;
            s_tot = s_Curve + 2 * s_Straight;

            size_chi = size(chi, 1);
            x = eta * 0;
            y = eta * 0;

            for i = 1:size_chi
                if (0 <= chi_arc(i,1)) && (chi_arc(i,1) < s_Straight/s_tot)
                    % Right straight segment
                    lambda = s_Straight/s_tot - chi_arc(i,1);
                    x_body = eta_phys(i,:)*sin(s.curvilinear_mapping.theta) ...
                        + s.curvilinear_mapping.R*sin(s.curvilinear_mapping.theta) ...
                        + lambda*cos(s.curvilinear_mapping.theta)*s_tot;
                    y_body = eta_phys(i,:)*cos(s.curvilinear_mapping.theta) ...
                        + (s.curvilinear_mapping.L - s.curvilinear_mapping.R) ...
                        + s.curvilinear_mapping.R*cos(s.curvilinear_mapping.theta) ...
                        - lambda*sin(s.curvilinear_mapping.theta)*s_tot;
                elseif (s_Straight/s_tot <= chi_arc(i,1)) && (chi_arc(i,1) < (s_Straight+s_Curve)/s_tot)
                    % Circular nose cap
                    angle = (chi_arc(i,1) - s_Straight/s_tot) * 2*s.curvilinear_mapping.theta*s_tot/s_Curve ...
                        + pi/2 - s.curvilinear_mapping.theta;
                    x_body = (eta_phys(i,:) + s.curvilinear_mapping.R) .* cos(angle);
                    y_body = (eta_phys(i,:) + s.curvilinear_mapping.R) .* sin(angle) ...
                        - s.curvilinear_mapping.R + s.curvilinear_mapping.L;
                else
                    % Left straight segment
                    lambda = chi_arc(i,1) - (s_Straight + s_Curve)/s_tot;
                    x_body = -eta_phys(i,:)*sin(s.curvilinear_mapping.theta) ...
                        - s.curvilinear_mapping.R*sin(s.curvilinear_mapping.theta) ...
                        - lambda*cos(s.curvilinear_mapping.theta)*s_tot;
                    y_body =  eta_phys(i,:)*cos(s.curvilinear_mapping.theta) ...
                        + (s.curvilinear_mapping.L - s.curvilinear_mapping.R) ...
                        + s.curvilinear_mapping.R*cos(s.curvilinear_mapping.theta) ...
                        - lambda*sin(s.curvilinear_mapping.theta)*s_tot;
                end
                % Body-aligned to physical frame: x = -y_body, y = x_body
                x(i,:) = -y_body;
                y(i,:) =  x_body;
            end

        %% Channel (with optional channel angle)
        case "channel"
            % Wall-normal refinement: cluster eta points near walls (0 and 1)
            % p = eta_refinement_power controls intensity (1 = none, 2 = quadratic, >2 = stronger)
            p = s.curvilinear_mapping.eta_refinement_power;
            eta_r = zeros(size(eta));
            mask = eta <= 0.5;
            eta_r(mask)  = 2^(p-1) * eta(mask).^p;
            eta_r(~mask) = 1 - 2^(p-1) * (1 - eta(~mask)).^p;

            x = chi .* cos(s.curvilinear_mapping.channel_angle) * s.curvilinear_mapping.Lx - eta_r .* sin(s.curvilinear_mapping.channel_angle) * s.curvilinear_mapping.Ly;
            y = eta_r .* cos(s.curvilinear_mapping.channel_angle) * s.curvilinear_mapping.Ly + chi .* sin(s.curvilinear_mapping.channel_angle) * s.curvilinear_mapping.Lx;

        %% Lid-driven cavity (nonlinear LDC refinement mapping)
        case "lid_driven_cavity"
            r = s.curvilinear_mapping.refinement_ratio_LDC;
            f_max = 1 + 2*(r-1)/3;

            chi_phys = (chi + 4*(r-1)*(chi.^2/2 - chi.^3/3)) / f_max;
            eta_phys = (eta + 4*(r-1)*(eta.^2/2 - eta.^3/3)) / f_max;

            x = -eta_phys;
            y = 1 - chi_phys;

        %% Wedge / Linear Interaction Analysis (linear distance falloff)
        case "wedge"
            distance_max = s.curvilinear_mapping.dRs ...
                + (s.curvilinear_mapping.dRe - s.curvilinear_mapping.dRs) * (1 - chi);
            x = -eta .* distance_max;
            y = s.curvilinear_mapping.Lx * (1 - chi);

        otherwise
            error('GO_TO_PHYSICAL_SPACE type not defined')
    end
end
