function s = INITIAL_SOLUTION(s, chemistry)
% INITIAL_SOLUTION  Set initial flow field to uniform freestream conditions.
%
%   s = INITIAL_SOLUTION(s, chemistry)
%
%   Initializes all conservative flow variables (density, momentum, and
%   total energy) to the upstream freestream values across the entire
%   computational domain. After setting the uniform field, the chemistry
%   equilibrium state is updated. For frozen (non-equilibrium) chemistry,
%   gamma_star and cv_star are also evaluated from the thermodynamic tables.
%
%   Inputs:
%       s  - (struct) Solver data structure. Must contain at minimum:
%           .curvilinear_mapping.boundary_type  - (string) Geometry type:
%                                                 'channel', 'MSL', 'blunt_cone',
%                                                 'circle', 'lid_driven_cavity'
%           .var.rho_E                  - (matrix) Total energy field (for sizing)
%           .freestream.rho_u_0         - (scalar) Non-dim upstream x-momentum
%           .freestream.rho_v_0         - (scalar) Non-dim upstream y-momentum
%           .freestream.rho_0           - (scalar) Non-dim upstream density
%           .freestream.rho_E_0         - (scalar) Non-dim upstream total energy
%           .freestream.rho_factor      - (scalar) Density scaling factor
%           .freestream.energy_factor   - (scalar) Energy scaling factor
%           .freestream.cv              - (scalar) Freestream cv for normalization
%           .chemistry.is_chemistry_enabled  - (logical) Chemistry toggle
%           .chemistry.chemical_equilibrium  - (logical) Equilibrium model flag
%           .mesh.y_Ext                 - (matrix) Extended y-coordinates (channel)
%           .curvilinear_mapping.Ly     - (scalar) Channel height (channel case)
%       chemistry - (struct) Chemistry model with evaluation functions:
%                     .frozen.eval_gamma_star - effective gamma evaluator
%                     .frozen.eval_cv_star    - effective cv evaluator
%
%   Outputs:
%       s  - (struct) Updated solver structure with uniform initial
%                   flow field and chemistry properties.
%
%   Notes:
%       - All supported boundary types receive the same uniform
%         initialization; the branching exists for future extensibility.
%       - The multiplication by zero (rho_E*0) is used to create an array
%         of the correct size filled with the upstream value.
%
%   See also: SET_FREESTREAM_PROPERTIES, UPDATE_CHEMISTRY_EQUILIBRIUM
%
%   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Set freestream reference properties
    s = SET_FREESTREAM_PROPERTIES(s, chemistry);

    %% Initialize conservative variables to uniform upstream state
    field_size = size(s.var.rho_E);

    if strcmp(s.curvilinear_mapping.boundary_type, 'MSL')           || ...
       strcmp(s.curvilinear_mapping.boundary_type, 'blunt_cone')    || ...
       strcmp(s.curvilinear_mapping.boundary_type, 'circle')        || ...
       strcmp(s.curvilinear_mapping.boundary_type, 'lid_driven_cavity')
        s.var.rho_u = zeros(field_size) + s.freestream.rho_u_0;
        s.var.rho_v = zeros(field_size) + s.freestream.rho_v_0;
        s.var.rho   = zeros(field_size) + s.freestream.rho_0;
        s.var.rho_E = zeros(field_size) + s.freestream.rho_E_0;
    elseif strcmp(s.curvilinear_mapping.boundary_type, 'channel')
        % Parabolic velocity profile: u(y) = u_max * y*(2-y)
        %   zero at y = 0 and y = 2, maximum at y = 1
        y = 2 * s.mesh.y_Ext ./ s.curvilinear_mapping.Ly;  % Normalize to [0,2] for the profile
        u_profile = y .* (2 - y);  % ranges from 0 to 1

        s.var.rho   = zeros(field_size) + s.freestream.rho_0;
        s.var.rho_v = zeros(field_size);  % v = 0
        s.var.rho_u = s.freestream.rho_u_0 * u_profile;

        % Internal energy from freestream, then add local kinetic energy
        rho_e_internal = s.freestream.rho_E_0 ...
            - 0.5 * (s.freestream.rho_u_0^2 + s.freestream.rho_v_0^2) / s.freestream.rho_0;
        s.var.rho_E = rho_e_internal + 0.5 * s.var.rho_u.^2 ./ s.var.rho;
    else
        error('Unsupported boundary type: %s', s.curvilinear_mapping.boundary_type);
    end

    %% Update chemistry equilibrium state
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);

    %% Evaluate frozen chemistry properties if applicable
    if s.chemistry.is_chemistry_enabled
        if ~s.chemistry.chemical_equilibrium
            rho = s.var.rho;
            e = (s.var.rho_E - (s.var.rho_u.^2 + s.var.rho_v.^2) ./ rho / 2) ./ rho;
            s.var.gamma_star = chemistry.frozen.eval_gamma_star( ...
                s.freestream.rho_factor * rho, ...
                s.freestream.energy_factor * e);
            s.var.cv_star = chemistry.frozen.eval_cv_star( ...
                s.freestream.rho_factor * rho, ...
                s.freestream.energy_factor * e) ./ s.freestream.cv;
        end
    end
end
