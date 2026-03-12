% input_file.m - Input configuration for 3D axisymmetric blunted cone (MSL-type)
%
% Description:
%   Defines all simulation parameters for a hypersonic shock instability
%   analysis over a blunted cone (MSL-type) geometry. This configuration
%   file is loaded by main.m via LOAD_INPUT_VARIABLES and populates the
%   'solution' struct.
%
% Usage:
%   This file is evaluated by LOAD_INPUT_VARIABLES and should not be called
%   directly. All parameters are written into the 'solution' struct.
%
% Notes:
%   - When is_chemistry_enabled is false, the solver uses gamma/Mach/Re/Pr.
%     When true, freestream conditions (u, rho, T, etc.) are used instead.
%
% Part of: Hypersonics Stability MATLAB Solver - Test Cases Module

%% ========================================================================
%  Solver Options
%  ========================================================================

solution.restart               = false;              % Resume from a previous run
solution.remesh                = true;               % Regenerate the mesh on restart
solution.restart_from_file     = false;              % Load restart data from a .mat file
solution.filename_restart         = "./solution.mat";   % Path to .mat restart file

if ~exist('solution_save', 'var') || isempty(solution_save)
    solution_save = [];
end

%% ========================================================================
%  Simulation Parameters
%  ========================================================================

%% PDE dimension
solution.PDE_dimension = "3D-axisymmetric";  % Options: "2D", "3D-axisymmetric"

%% Chemistry model
solution.chemistry.is_chemistry_enabled       = true;             % Enable thermochemistry (false = calorically perfect gas)
solution.chemistry.chemistry_type        = "Chemical-RTVE";  % Model: "Frozen-RTV", "Frozen-RTVE", "Chemical-RTV", "Chemical-RTVE"
solution.chemistry.chemical_equilibrium  = true;             % Use equilibrium chemistry
solution.chemistry.non_equilibrium_model = "linear";         % Non-equilibrium model: "linear", "quadratic"
solution.chemistry.chemistry_composition      = "Earth";          % Atmospheric composition: "Earth", "Mars", "CO2"

%% Mesh parameters
solution.mesh.Nchi = 400;  % Grid points along chi (wall boundary)
solution.mesh.Neta = 40;   % Grid points along eta (wall-normal)

%% Geometry: blunt cone (MSL-type)
solution.curvilinear_mapping.boundary_type         = "blunt_cone";  % Geometry type
solution.curvilinear_mapping.theta                 = 70 * pi / 180; % Half-angle of cone tip [rad]
solution.curvilinear_mapping.R                     = 1;             % Tip sphere radius
solution.curvilinear_mapping.L                     = 3;             % x distance from outflow to tip
solution.curvilinear_mapping.dRe                   = 0.6;           % Wall distance of domain boundary at outflow edge
solution.curvilinear_mapping.dRs                   = 0.4;           % Wall distance of domain boundary at stagnation
solution.curvilinear_mapping.refinement_stagnation.state        = true;   % Stagnation region mesh refinement
solution.curvilinear_mapping.refinement_stagnation.BL_thickness = 0.2;    % Thickness of refined region
solution.curvilinear_mapping.refinement_stagnation.intensity    = 0.99;   % Refinement intensity (0 = none, 1 = max)
solution.curvilinear_mapping.refinement_wall.state              = false;  % Wall region mesh refinement
solution.curvilinear_mapping.smooth_mesh = true;  % Smooth mesh for curvature discontinuities

%% Boundary conditions
%         _____eta1_______
%        |                |
%   chi0 |                | chi1
%        |                |
%        ------eta0--------
%
% Available types:
%   'inflow_subsonic'  'inflow_supersonic'  'periodic'  'shock'
%   'outflow_supersonic'  'outflow_subsonic'  'outflow_NRCBC'
%   'no_slip_adiabatic'  'no_slip_isothermal'  'symmetry'

solution.boundary_conditions.boundary_eta0.name        = 'no_slip_adiabatic';  % Wall BC
solution.boundary_conditions.boundary_eta1.name        = 'shock';              % Upper (shock) BC
solution.boundary_conditions.boundary_chi0.name        = 'symmetry';           % Inflow/symmetry BC
solution.boundary_conditions.boundary_chi1.name        = 'outflow_NRCBC';      % Outflow BC

%% Freestream conditions (chemistry mode)
% Earth atmosphere at 50 km
solution.freestream.u   = 4000;      % Freestream velocity [m/s]
solution.freestream.v   = 0;         % Transverse velocity [m/s]
solution.freestream.rho = 0.001;     % Freestream density [kg/m^3]
solution.freestream.T   = 270;       % Freestream temperature [K]
solution.freestream.Re  = 100000;    % Reynolds number

%% Time integration
solution.time_integration.N_iter            = 1000;           % Number of time steps
solution.time_integration.time_integrator   = "Explicit_RK4"; % Scheme: "Explicit_RK4", "Implicit_Euler"
solution.time_integration.CFL               = 2;             % CFL number for adaptive time stepping
solution.time_integration.dt                = 0.000001;      % Initial time step [s]
solution.time_integration.max_dt            = 0.3;           % Maximum allowable time step [s]

%% Shock fitting
solution.shock.enabled                  = true;         % Enable shock fitting
solution.shock.feedback                 = true;         % Dynamic shock response to perturbations
solution.shock.interpolate              = "2nd";        % Extrapolation order: "1st", "2nd", "3rd"
solution.shock.initial_shock_dist       = 0.3;          % Initial shock standoff distance from wall
solution.shock.remesh_shock_distance    = 1.2;          % Domain boundary distance relative to shock position
solution.shock.relaxation               = 1.0;          % Shock speed relaxation (0 = frozen, 1 = physical, <1 = under-relaxed)
solution.shock.formulation              = "Lagrangian"; % Frame formulation: "Lagrangian", "Eulerian"
solution.shock.fitting                  = "csaps";      % Spline method: "csaps" (cubic smoothing spline)
solution.shock.spline_param             = 1 - 1e-8 * 100^3 / solution.mesh.Nchi^3; % Smoothing parameter for csaps
solution.shock.initial_beta             = 55 * pi / 180; % Initial guess for shock angle initialization [rad]

%% Stability analysis and transient growth
solution.stability_analysis.perturbation_magnitude  = 1e-8;                    % Linearization perturbation size
solution.stability_analysis.eigenvalue_solver       = "GPU_TIMESTEPPER_ARNOLDI"; % Eigensolver: "CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"

%% Live plots during time-marching
solution.running_plot.enabled   = true;                       % Show live plots during simulation
solution.running_plot.variable  = {'div_U', 'rho', 'vort'};  % Variables to plot: 'rho', 'u', 'v', 'p', 'T', 'div_U', 'vort', etc.
solution.running_plot.timesteps = 100;                       % Update interval [time steps]

fprintf("Input data processed. \n")
