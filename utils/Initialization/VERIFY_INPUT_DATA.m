

function s = VERIFY_INPUT_DATA(s)
% VERIFY_INPUT_DATA  Validate and augment the solver input data structure.
%
%   s = VERIFY_INPUT_DATA(s)
%
%   Performs comprehensive validation of all solver configuration fields,
%   checks types and allowed values, enforces consistency constraints
%   (e.g., periodic BCs must appear on both sides), and sets default values
%   for optional parameters that are not explicitly provided.
%
%   Inputs:
%       s - (struct) Solver configuration structure with the following
%           groups of fields (see Notes for details):
%           .restart, .remesh              - (logical) Solver control flags
%           .PDE_dimension                 - (string) "2D" or "3D-axisymmetric"
%           .chemistry.is_chemistry_enabled      - (logical) Chemistry toggle
%           .chemistry.chemistry_type            - (string) Chemistry model
%           .chemistry.chemical_equilibrium      - (logical) Equilibrium flag
%           .chemistry.non_equilibrium_model     - (string) NEQ model type
%           .chemistry.chemistry_composition     - (string) Planet atmosphere
%           .freestream.u, .freestream.v         - (double) Velocity components
%           .freestream.rho, .freestream.T       - (double) For chemistry mode
%           .freestream.gamma, .freestream.Mach  - (double) For ideal-gas mode
%           .freestream.Re or .freestream.L_ref  - (double) Length scale
%           .mesh.Nchi, .mesh.Neta               - (int) Grid dimensions
%           .curvilinear_mapping.boundary_type    - (string) Geometry type
%           .boundary_conditions.boundary_*      - (struct) BC specifications
%           .time_integration.*                  - Time integrator settings
%           .numerical_dissipation.*             - (optional) Dissipation coeffs
%           .shock.*                             - (optional) Shock-fitting params
%           .stability_analysis.*                - (optional) Stability settings
%           .running_plot.*                      - (optional) Live plot settings
%
%   Outputs:
%       s - (struct) Validated structure with defaults set for optional
%           fields (numerical_dissipation, refinement, freestream
%           disturbance, stability_analysis.perturb_shock, etc.)
%
%   Notes:
%       - Chemistry-disabled mode sets defaults: chemistry_type = "None",
%         chemical_equilibrium = false, chemistry_composition = "None".
%       - Periodic BCs require matching on opposite boundaries and grid
%         dimensions divisible by 3.
%       - Supported boundary types: "circle", "lid_driven_cavity",
%         "channel", "MSL", "blunt_cone".
%       - Supported BC names: 'inflow_subsonic', 'inflow_supersonic',
%         'periodic', 'shock', 'outflow_supersonic', 'outflow_subsonic',
%         'outflow_NRCBC', 'no_slip_adiabatic', 'no_slip_isothermal',
%         'symmetry'.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    % 1. Solver Options
    assert(isfield(s, 'restart'), 'Missing field: restart');
    assert(islogical(s.restart), 'restart must be a logical');
    assert(isfield(s, 'remesh'), 'Missing field: remesh');
    assert(islogical(s.remesh), 'remesh must be a logical');
    if s.restart
        assert(isfield(s, 'restart_from_file'), 'Missing field: restart_from_file');
        assert(islogical(s.restart_from_file), 'restart_from_file must be a logical');
        if s.restart_from_file
            assert(isfield(s, 'filename_restart'), 'Missing field: filename_restart');
            assert(isstring(s.filename_restart) || ischar(s.filename_restart), 'filename_restart must be a string or char');
        end
    end

    % 2. PDE dimension
    assert(isfield(s, 'PDE_dimension'), 'Missing field: PDE_dimension');
    assert(ismember(s.PDE_dimension, ["2D", "3D-axisymmetric"]), 'PDE_dimension must be "2D" or "3D-axisymmetric"');

    % 3. Chemistry
    assert(isfield(s, 'chemistry'), 'Missing field: chemistry');
    assert(isfield(s.chemistry, 'is_chemistry_enabled'), 'Missing field: chemistry.is_chemistry_enabled');
    assert(islogical(s.chemistry.is_chemistry_enabled), 'chemistry.is_chemistry_enabled must be a logical');
    if s.chemistry.is_chemistry_enabled
        assert(isfield(s.chemistry, 'chemistry_type'), 'Missing field: chemistry.chemistry_type');
        assert(ismember(s.chemistry.chemistry_type, ["Frozen-RTV", "Frozen-RTVE", "Chemical-RTV", "Chemical-RTVE"]), 'Invalid chemistry_type');
        assert(isfield(s.chemistry, 'chemical_equilibrium'), 'Missing field: chemistry.chemical_equilibrium');
        assert(islogical(s.chemistry.chemical_equilibrium), 'chemistry.chemical_equilibrium must be a logical');
        assert(isfield(s.chemistry, 'non_equilibrium_model'), 'Missing field: chemistry.non_equilibrium_model');
        assert(ismember(s.chemistry.non_equilibrium_model, ["linear", "quadratic"]), 'Invalid non_equilibrium_model');
        assert(isfield(s.chemistry, 'chemistry_composition'), 'Missing field: chemistry.chemistry_composition');
        assert(ismember(s.chemistry.chemistry_composition, ["Earth", "Mars", "CO2"]), 'Invalid chemistry_composition');
    else
        s.chemistry.is_chemistry_enabled  = false;
        s.chemistry.chemistry_type   = "None";
        s.chemistry.chemical_equilibrium = false;
        s.chemistry.chemistry_composition = "None";
        s.chemistry.planet = "NONE";
        s.chemistry.gas_model = "NONE";
    end

    % 4. Freestream
    assert(isfield(s, 'freestream'), 'Missing field: freestream');
    assert(isfield(s.freestream, 'u'), 'Missing field: freestream.u');
    assert(isnumeric(s.freestream.u), 'freestream.u must be numeric');
    assert(isfield(s.freestream, 'v'), 'Missing field: freestream.v');
    assert(isnumeric(s.freestream.v), 'freestream.v must be numeric');
    if s.chemistry.is_chemistry_enabled
        assert(isfield(s.freestream, 'rho'), 'Missing field: freestream.rho');
        assert(isnumeric(s.freestream.rho) && s.freestream.rho > 0, 'freestream.rho must be positive numeric');
        assert(isfield(s.freestream, 'T'), 'Missing field: freestream.T');
        assert(isnumeric(s.freestream.T) && s.freestream.T > 0, 'freestream.T must be positive numeric');
    else
        assert(isfield(s.freestream, 'gamma'), 'Missing field: freestream.gamma');
        assert(isnumeric(s.freestream.gamma) && s.freestream.gamma > 1, 'freestream.gamma must be > 1');
        assert(isfield(s.freestream, 'Mach'), 'Missing field: freestream.Mach');
        assert(isnumeric(s.freestream.Mach) && s.freestream.Mach > 0, 'freestream.Mach must be positive numeric');
        assert(isfield(s.freestream, 'Pr'), 'Missing field: freestream.Pr');
        assert(isnumeric(s.freestream.Pr) && s.freestream.Pr > 0, 'freestream.Pr must be positive numeric');
    end
    
    if s.chemistry.is_chemistry_enabled
        has_L = isfield(s.freestream, 'L_ref');
        has_Re = isfield(s.freestream, 'Re');
        assert(xor(has_L, has_Re), 'Exactly one of freestream.L_ref or freestream.Re must be specified, not both.');
        
        if has_L
            assert(isnumeric(s.freestream.L_ref) && s.freestream.L_ref > 0, 'freestream.L_ref must be positive numeric');
        end
        if has_Re
            assert(isnumeric(s.freestream.Re) && s.freestream.Re > 0, 'freestream.Re must be positive numeric');
        end
    else
        assert(isfield(s.freestream, 'Re'), 'Missing field: freestream.Re');
        assert(isnumeric(s.freestream.Re) && s.freestream.Re > 0, 'freestream.Re must be positive numeric');
        assert(~isfield(s.freestream, 'L_ref'), 'freestream.L_ref should not be specified when chemistry is disabled. Specify freestream.Re instead.');
        if ~isfield(s.freestream, 'rho')
            s.freestream.rho = 1;
        end
        if ~isfield(s.freestream, 'T')
            s.freestream.T = 1;
        end
    end

    % 4b. Freestream Perturbation (optional - initialize to defaults if absent)
    if ~isfield(s.freestream, 'disturbance') || isempty(s.freestream.disturbance)
        s.freestream.disturbance.k_x       = 1;
        s.freestream.disturbance.k_y       = 1;
        s.freestream.disturbance.amplitude = [0, 0, 0, 0];
    else
        if ~isfield(s.freestream.disturbance, 'k_x')
            s.freestream.disturbance.k_x = 1;
        else
            assert(isnumeric(s.freestream.disturbance.k_x), 'freestream.disturbance.k_x must be numeric');
        end
        if ~isfield(s.freestream.disturbance, 'k_y')
            s.freestream.disturbance.k_y = 1;
        else
            assert(isnumeric(s.freestream.disturbance.k_y), 'freestream.disturbance.k_y must be numeric');
        end
        if ~isfield(s.freestream.disturbance, 'amplitude')
            s.freestream.disturbance.amplitude = [0, 0, 0, 0];
        else
            assert(isnumeric(s.freestream.disturbance.amplitude) && numel(s.freestream.disturbance.amplitude) == 4, ...
                'freestream.disturbance.amplitude must be a numeric vector of length 4 [rho, u, v, rhoE]');
        end
    end

    % 5. Mesh
    assert(isfield(s, 'mesh'), 'Missing field: mesh');
    assert(isfield(s.mesh, 'Nchi'), 'Missing field: mesh.Nchi');
    assert(isnumeric(s.mesh.Nchi) && s.mesh.Nchi > 0 && mod(s.mesh.Nchi, 1) == 0, 'mesh.Nchi must be a positive integer');
    assert(isfield(s.mesh, 'Neta'), 'Missing field: mesh.Neta');
    assert(isnumeric(s.mesh.Neta) && s.mesh.Neta > 0 && mod(s.mesh.Neta, 1) == 0, 'mesh.Neta must be a positive integer');

    % 6. Curvilinear Mapping
    assert(isfield(s, 'curvilinear_mapping'), 'Missing field: curvilinear_mapping');
    assert(isfield(s.curvilinear_mapping, 'boundary_type'), 'Missing field: curvilinear_mapping.boundary_type');
    assert(ismember(s.curvilinear_mapping.boundary_type, ["circle", "lid_driven_cavity", "channel", "MSL", "blunt_cone"]), 'Invalid boundary_type');

    if ~isfield(s.curvilinear_mapping, 'refinement_stagnation') || ~isfield(s.curvilinear_mapping.refinement_stagnation, 'state')
        s.curvilinear_mapping.refinement_stagnation.state = false;
    end
    if ~isfield(s.curvilinear_mapping, 'refinement_wall') || ~isfield(s.curvilinear_mapping.refinement_wall, 'state')
        s.curvilinear_mapping.refinement_wall.state = false;
    end

    % 7. Boundary Conditions
    assert(isfield(s, 'boundary_conditions'), 'Missing field: boundary_conditions');
    valid_bcs = {'inflow_subsonic', 'inflow_supersonic', 'periodic', 'shock', 'outflow_supersonic', 'outflow_subsonic', 'outflow_NRCBC', 'no_slip_adiabatic', 'no_slip_isothermal', 'symmetry'};
    assert(isfield(s.boundary_conditions, 'boundary_eta0'), 'Missing field: boundary_conditions.boundary_eta0.name');
    assert(ismember(s.boundary_conditions.boundary_eta0.name, valid_bcs), 'Invalid boundary_eta0');
    assert(isfield(s.boundary_conditions, 'boundary_eta1'), 'Missing field: boundary_conditions.boundary_eta1.name');
    assert(ismember(s.boundary_conditions.boundary_eta1.name, valid_bcs), 'Invalid boundary_eta1');
    assert(isfield(s.boundary_conditions, 'boundary_chi0'), 'Missing field: boundary_conditions.boundary_chi0.name');
    assert(ismember(s.boundary_conditions.boundary_chi0.name, valid_bcs), 'Invalid boundary_chi0');
    assert(isfield(s.boundary_conditions, 'boundary_chi1'), 'Missing field: boundary_conditions.boundary_chi1.name');
    assert(ismember(s.boundary_conditions.boundary_chi1.name, valid_bcs), 'Invalid boundary_chi1');

    % 7b. Periodic boundary condition consistency
    %     If one side of a direction is periodic, the opposite side must also be periodic.
    eta0_periodic = strcmp(s.boundary_conditions.boundary_eta0.name, 'periodic');
    eta1_periodic = strcmp(s.boundary_conditions.boundary_eta1.name, 'periodic');
    assert(eta0_periodic == eta1_periodic, ...
        'Periodic boundary conditions must be set on both sides: boundary_eta0 and boundary_eta1 must both be periodic or both non-periodic.');

    chi0_periodic = strcmp(s.boundary_conditions.boundary_chi0.name, 'periodic');
    chi1_periodic = strcmp(s.boundary_conditions.boundary_chi1.name, 'periodic');
    assert(chi0_periodic == chi1_periodic, ...
        'Periodic boundary conditions must be set on both sides: boundary_chi0 and boundary_chi1 must both be periodic or both non-periodic.');

    % 7c. When both eta boundaries are periodic, Neta must be a multiple of 3
    if eta0_periodic && eta1_periodic
        assert(mod(s.mesh.Neta, 3) == 0, ...
            'When both boundary_eta0 and boundary_eta1 are periodic, mesh.Neta must be a multiple of 3.');
    end
    if chi0_periodic && chi1_periodic
        assert(mod(s.mesh.Nchi, 3) == 0, ...
            'When both boundary_chi0 and boundary_chi1 are periodic, mesh.Nchi must be a multiple of 3.');
    end

    % 8. Time Integration
    assert(isfield(s, 'time_integration'), 'Missing field: time_integration');
    assert(isfield(s.time_integration, 'N_iter'), 'Missing field: time_integration.N_iter');
    assert(isnumeric(s.time_integration.N_iter) && s.time_integration.N_iter > 0 && mod(s.time_integration.N_iter, 1) == 0, 'time_integration.N_iter must be a positive integer');
    assert(isfield(s.time_integration, 'time_integrator'), 'Missing field: time_integration.time_integrator');
    assert(ismember(s.time_integration.time_integrator, ["Explicit_RK4", "Implicit_Euler"]), 'Invalid time_integrator');
    if s.time_integration.time_integrator == "Implicit_Euler"
        assert(isfield(s.time_integration, 'tolerance'), 'Missing field: time_integration.tolerance for Implicit_Euler');
        assert(isnumeric(s.time_integration.tolerance) && s.time_integration.tolerance > 0, 'time_integration.tolerance must be positive numeric');
        assert(isfield(s.time_integration, 'max_iter_implicit'), 'Missing field: time_integration.max_iter_implicit for Implicit_Euler');
        assert(isnumeric(s.time_integration.max_iter_implicit) && s.time_integration.max_iter_implicit > 0 && mod(s.time_integration.max_iter_implicit, 1) == 0, 'time_integration.max_iter_implicit must be a positive integer');
        assert(isfield(s.time_integration, 'relax_factor'), 'Missing field: time_integration.relax_factor for Implicit_Euler');
        assert(isnumeric(s.time_integration.relax_factor) && s.time_integration.relax_factor > 0 && s.time_integration.relax_factor <= 1, 'time_integration.relax_factor must be between 0 and 1');
    end
    assert(isfield(s.time_integration, 'CFL'), 'Missing field: time_integration.CFL');
    assert(isnumeric(s.time_integration.CFL) && s.time_integration.CFL > 0, 'time_integration.CFL must be positive numeric');
    assert(isfield(s.time_integration, 'dt'), 'Missing field: time_integration.dt');
    assert(isnumeric(s.time_integration.dt) && s.time_integration.dt > 0, 'time_integration.dt must be positive numeric');
    assert(isfield(s.time_integration, 'max_dt'), 'Missing field: time_integration.max_dt');
    assert(isnumeric(s.time_integration.max_dt) && s.time_integration.max_dt > 0, 'time_integration.max_dt must be positive numeric');
    if ~isfield(s.time_integration, 'plot_residual')
        s.time_integration.plot_residual = false;
    else

    % 9. Numerical Dissipation
    if ~isfield(s, 'numerical_dissipation')
        s.numerical_dissipation.mu_rho = 0;
        s.numerical_dissipation.mu_rho_u = 0;
        s.numerical_dissipation.mu_rho_v = 0;
        s.numerical_dissipation.mu_rho_E = 0;
        s.numerical_dissipation.rho_min = 0.05;
    else
        assert(isfield(s.numerical_dissipation, 'mu_rho'), 'Missing field: numerical_dissipation.mu_rho');
        assert(isnumeric(s.numerical_dissipation.mu_rho) && s.numerical_dissipation.mu_rho >= 0, 'numerical_dissipation.mu_rho must be non-negative numeric');
        assert(isfield(s.numerical_dissipation, 'mu_rho_u'), 'Missing field: numerical_dissipation.mu_rho_u');
        assert(isnumeric(s.numerical_dissipation.mu_rho_u) && s.numerical_dissipation.mu_rho_u >= 0, 'numerical_dissipation.mu_rho_u must be non-negative numeric');
        assert(isfield(s.numerical_dissipation, 'mu_rho_v'), 'Missing field: numerical_dissipation.mu_rho_v');
        assert(isnumeric(s.numerical_dissipation.mu_rho_v) && s.numerical_dissipation.mu_rho_v >= 0, 'numerical_dissipation.mu_rho_v must be non-negative numeric');
        assert(isfield(s.numerical_dissipation, 'mu_rho_E'), 'Missing field: numerical_dissipation.mu_rho_E');
        assert(isnumeric(s.numerical_dissipation.mu_rho_E) && s.numerical_dissipation.mu_rho_E >= 0, 'numerical_dissipation.mu_rho_E must be non-negative numeric');
        assert(isfield(s.numerical_dissipation, 'rho_min'), 'Missing field: numerical_dissipation.rho_min');
        assert(isnumeric(s.numerical_dissipation.rho_min) && s.numerical_dissipation.rho_min >= 0, 'numerical_dissipation.rho_min must be non-negative numeric');
    end

    % 10. Shock Fitting
    assert(isfield(s, 'shock'), 'Missing field: shock');
    assert(isfield(s.shock, 'enabled'), 'Missing field: shock.enabled');
    assert(islogical(s.shock.enabled), 'shock.enabled must be a logical');
    if s.shock.enabled
        assert(isfield(s.shock, 'feedback'), 'Missing field: shock.feedback');
        assert(islogical(s.shock.feedback), 'shock.feedback must be a logical');
        assert(isfield(s.shock, 'interpolate'), 'Missing field: shock.interpolate');
        assert(ismember(s.shock.interpolate, ["1st", "2nd", "3rd"]), 'Invalid shock.interpolate');
        assert(isfield(s.shock, 'initial_shock_dist'), 'Missing field: shock.initial_shock_dist');
        assert(isnumeric(s.shock.initial_shock_dist) && s.shock.initial_shock_dist > 0, 'shock.initial_shock_dist must be positive numeric');
        assert(isfield(s.shock, 'remesh_shock_distance'), 'Missing field: shock.remesh_shock_distance');
        assert(isnumeric(s.shock.remesh_shock_distance) && s.shock.remesh_shock_distance > 0, 'shock.remesh_shock_distance must be positive numeric');
        assert(isfield(s.shock, 'relaxation'), 'Missing field: shock.relaxation');
        assert(isnumeric(s.shock.relaxation) && s.shock.relaxation >= 0 && s.shock.relaxation <= 1, 'shock.relaxation must be between 0 and 1');
        assert(isfield(s.shock, 'formulation'), 'Missing field: shock.formulation');
        assert(ismember(s.shock.formulation, ["Lagrangian", "Eulerian"]), 'Invalid shock.formulation');
        assert(isfield(s.shock, 'fitting'), 'Missing field: shock.fitting');
        assert(ismember(s.shock.fitting, ["csaps"]), 'Invalid shock.fitting: only "csaps" is supported');
        assert(isfield(s.shock, 'spline_param'), 'Missing field: shock.spline_param');
        assert(isnumeric(s.shock.spline_param) && s.shock.spline_param >= 0 && s.shock.spline_param <= 1, 'shock.spline_param must be between 0 and 1');
    end

    % 11. Stability Analysis
    assert(isfield(s, 'stability_analysis'), 'Missing field: stability_analysis');
    assert(isfield(s.stability_analysis, 'perturbation_magnitude'), 'Missing field: stability_analysis.perturbation_magnitude');
    assert(isnumeric(s.stability_analysis.perturbation_magnitude) && s.stability_analysis.perturbation_magnitude > 0, 'stability_analysis.perturbation_magnitude must be positive numeric');
    assert(isfield(s.stability_analysis, 'eigenvalue_solver'), 'Missing field: stability_analysis.eigenvalue_solver');
    assert(ismember(s.stability_analysis.eigenvalue_solver, ["CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"]), 'Invalid eigenvalue_solver');
    if ~isfield(s.stability_analysis, 'perturb_shock')
        s.stability_analysis.perturb_shock = false;
    else
        assert(islogical(s.stability_analysis.perturb_shock), 'stability_analysis.perturb_shock must be a logical');
    end

    % 12. Running Plot
    assert(isfield(s, 'running_plot'), 'Missing field: running_plot');
    assert(isfield(s.running_plot, 'enabled'), 'Missing field: running_plot.enabled');
    assert(islogical(s.running_plot.enabled), 'running_plot.enabled must be a logical');
    if s.running_plot.enabled
        assert(isfield(s.running_plot, 'variable'), 'Missing field: running_plot.variable');
        assert(iscell(s.running_plot.variable), 'running_plot.variable must be a cell array of strings');
        assert(isfield(s.running_plot, 'timesteps'), 'Missing field: running_plot.timesteps');
        assert(isnumeric(s.running_plot.timesteps) && s.running_plot.timesteps > 0 && mod(s.running_plot.timesteps, 1) == 0, 'running_plot.timesteps must be a positive integer');
    end
end