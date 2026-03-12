function s = SHOCK_LINE_INITIALIZATION(s)
% SHOCK_LINE_INITIALIZATION  Initialize the shock line position for various geometries.
%   Sets the initial shock stand-off distance and shock point coordinates
%   in both physical (x, y) and element (chi, eta) space, depending on the
%   boundary type (circle, MSL, blunt_cone, or LIA).
%
%   s = SHOCK_LINE_INITIALIZATION(s)
%
%   Inputs:
%       s - Solution struct containing:
%           .restart                            - (logical) Skip if true
%           .mesh.Nchi                          - (int) Number of streamwise cells
%           .mesh.chi                           - (Nchi x Neta) element-space coords
%           .curvilinear_mapping.boundary_type  - (string) Geometry type
%           .curvilinear_mapping.R              - (double) Nose radius
%           .curvilinear_mapping.L              - (double) Body length (blunt_cone)
%           .curvilinear_mapping.dRs            - (double) Outer boundary at stagnation
%           .shock.initial_shock_dist           - (double) Shock standoff distance
%           .shock.initial_beta                 - (double) Initial shock angle (blunt_cone)
%
%   Outputs:
%       s - Updated struct with initialized shock point arrays:
%           .shock.points_x, .shock.points_y     - Physical coordinates
%           .shock.points_eta, .shock.points_chi  - Element-space coordinates
%
%   Notes:
%       - Only executes when s.restart is false (fresh start).
%       - Element-space coordinates (chi, eta) are in [0,1].
%       - For 'circle' and 'MSL', the shock eta is uniform at
%         initial_shock_dist / dRs.
%       - For 'blunt_cone', the shock follows the conical and curved
%         sections of the geometry with a constant stand-off distance.
%       - For 'wedge' (Linear Interaction Analysis), the shock is an
%         oblique line with angle shock_delta.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    if ~s.restart

        %% Common initialization
        s.shock.points_chi = s.mesh.chi(:, 1);
        s.shock.points_eta = 0 * s.mesh.chi(:, 1);

        if s.curvilinear_mapping.boundary_type == 'circle'
            %% Circular geometry: uniform eta in [0,1]
            s.shock.points_eta(:, 1) = s.shock.initial_shock_dist ...
                / s.curvilinear_mapping.dRs;
            [s.shock.points_x, s.shock.points_y] = GO_TO_PHYSICAL_SPACE( ...
                s.shock.points_chi, s.shock.points_eta, s);

        elseif s.curvilinear_mapping.boundary_type == 'MSL'
            %% MSL geometry: uniform eta in [0,1]
            s.shock.points_eta(:, 1) = s.shock.initial_shock_dist ...
                / s.curvilinear_mapping.dRs;
            [s.shock.points_x, s.shock.points_y] = GO_TO_PHYSICAL_SPACE( ...
                s.shock.points_chi, s.shock.points_eta, s);

        elseif s.curvilinear_mapping.boundary_type == 'blunt_cone'
            %% Blunt cone geometry: straight + curved sections
            s_Straight = ((s.curvilinear_mapping.L - s.curvilinear_mapping.R) + s.curvilinear_mapping.R * cos(s.shock.initial_beta)) ...
                / sin(s.shock.initial_beta);
            s_Curve = s.shock.initial_beta * s.curvilinear_mapping.R;
            s_tot = 2 * s_Curve + 2 * s_Straight;

            for i = 1:s.mesh.Nchi
                % Convert chi from [0,1] to arc-length coordinate (chi=0 left, chi=1 right)
                chi_arc_i = (1 - s.mesh.chi(i, 1)) / 2;

                if (chi_arc_i < s_Straight / s_tot)
                    % Straight conical section
                    s_len = ((s_Straight) / s_tot - chi_arc_i) * s_tot;
                    R = s.shock.initial_shock_dist + s.curvilinear_mapping.R;
                    x_body = R * sin(s.shock.initial_beta) ...
                        + cos(s.shock.initial_beta) * s_len;
                    y_body = R * cos(s.shock.initial_beta) ...
                        + (s.curvilinear_mapping.L - s.curvilinear_mapping.R) - sin(s.shock.initial_beta) * s_len;
                    % Body-aligned to physical frame
                    s.shock.points_x(i, 1) = -y_body;
                    s.shock.points_y(i, 1) = x_body;
                else
                    % Curved nose section
                    R = s.shock.initial_shock_dist + s.curvilinear_mapping.R;
                    lambda = ((s_Straight + s_Curve) / s_tot - chi_arc_i) / s.curvilinear_mapping.R;
                    x_body = R * sin(lambda);
                    y_body = R * cos(lambda) + (s.curvilinear_mapping.L - s.curvilinear_mapping.R);
                    % Body-aligned to physical frame
                    s.shock.points_x(i, 1) = -y_body;
                    s.shock.points_y(i, 1) = x_body;
                end
            end

            [s.shock.points_chi, s.shock.points_eta] = GO_TO_ELEMENT_SPACE( ...
                s.shock.points_x, s.shock.points_y, s);

        elseif s.curvilinear_mapping.boundary_type == 'wedge'
            %% Linear Interaction Analysis: oblique shock line
            if ~isfield("s", "shock_delta")
                s.shock_delta = pi / 8;
            end
            s.shock.points_y(:, 1) = s.mesh.y(:, 1);
            s.shock.points_x(:, 1) = -s.shock.chi0_shock ...
                - s.shock.points_y(:, 1) .* sin(s.shock_delta);
            [s.shock.points_chi, s.shock.points_eta] = GO_TO_ELEMENT_SPACE( ...
                s.shock.points_x, s.shock.points_y, s);
        end
    end
end
