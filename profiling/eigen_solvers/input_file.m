% input_file.m - Input configuration for hypersonic flow over a cylinder
%
% Description:
%   Defines all simulation parameters for a hypersonic shock instability
%   analysis over a circular cylinder geometry. This configuration file
%   is loaded by main.m via LOAD_INPUT_VARIABLES and populates the
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
solution.PDE_dimension = "2D";  % Options: "2D", "3D-axisymmetric"

%% Chemistry model
solution.chemistry.is_chemistry_enabled       = true;             % Enable thermochemistry (false = calorically perfect gas)
solution.chemistry.chemistry_type        = "Chemical-RTVE";  % Model: "Frozen-RTV", "Frozen-RTVE", "Chemical-RTV", "Chemical-RTVE"
solution.chemistry.chemical_equilibrium  = true;             % Use equilibrium chemistry
solution.chemistry.non_equilibrium_model = "linear";         % Non-equilibrium model: "linear", "quadratic"
solution.chemistry.chemistry_composition      = "Mars";           % Atmospheric composition: "Earth", "Mars", "CO2"

%% Mesh parameters
solution.mesh.Nchi = 200;  % Grid points along chi (wall boundary)
solution.mesh.Neta = 40;   % Grid points along eta (wall-normal)

%% Geometry: circle/cylinder
solution.curvilinear_mapping.boundary_type    = "circle";  % Geometry type
solution.curvilinear_mapping.R                = 1;         % Cylinder radius
solution.curvilinear_mapping.dRe              = 1.5;       % Wall distance of domain boundary at outflow edge
solution.curvilinear_mapping.dRs              = 0.5;       % Wall distance of domain boundary at stagnation
solution.curvilinear_mapping.circle_angle_extra = 0;       % Extra angular extent beyond 90 deg [rad]
solution.curvilinear_mapping.refinement_stagnation.state        = false;  % Stagnation region mesh refinement
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
solution.freestream.u   = 5000;      % Freestream velocity [m/s]
solution.freestream.v   = 0;         % Transverse velocity [m/s]
solution.freestream.rho = 0.001;     % Freestream density [kg/m^3]
solution.freestream.T   = 300;       % Freestream temperature [K]
solution.freestream.Re  = 10000;     % Reynolds number

%% Time integration
solution.time_integration.N_iter            = 1000;           % Number of time steps
solution.time_integration.time_integrator   = "Explicit_RK4"; % Scheme: "Explicit_RK4", "Implicit_Euler"
%solution.time_integration.tolerance         = 1e-3;          % Implicit solver convergence tolerance
%solution.time_integration.max_iter_implicit = 1000;          % Implicit solver max iterations
%solution.time_integration.relax_factor      = 0.9;           % Implicit solver under-relaxation factor
solution.time_integration.CFL               = 2;             % CFL number for adaptive time stepping
solution.time_integration.dt                = 0.000001;      % Initial time step [s]
solution.time_integration.max_dt            = 0.3;           % Maximum allowable time step [s]

%% Shock fitting
solution.shock.enabled                  = true;         % Enable shock fitting
solution.shock.feedback                 = true;         % Dynamic shock response to perturbations
solution.shock.interpolate              = "2nd";        % Extrapolation order: "1st", "2nd", "3rd"
solution.shock.initial_shock_dist       = 0.4;          % Initial shock standoff distance from wall
solution.shock.remesh_shock_distance    = 1.3;          % Domain boundary distance relative to shock position
solution.shock.relaxation               = 1.0;          % Shock speed relaxation (0 = frozen, 1 = physical, <1 = under-relaxed)
solution.shock.formulation              = "Lagrangian"; % Frame formulation: "Lagrangian", "Eulerian"
solution.shock.fitting                  = "csaps";      % Spline method: "csaps" (cubic smoothing spline)
solution.shock.spline_param             = 1 - 1e-6 * 100^3 / solution.mesh.Nchi^3; % Smoothing parameter for csaps

%% Stability analysis and transient growth
solution.stability_analysis.perturbation_magnitude  = 1e-8;                    % Linearization perturbation size
solution.stability_analysis.eigenvalue_solver       = "GPU_TIMESTEPPER_ARNOLDI"; % Eigensolver: "CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"

%% Live plots during time-marching
solution.running_plot.enabled   = true;                     % Show live plots during simulation
solution.running_plot.variable  = {'vort'};                 % Variables to plot: 'rho', 'u', 'v', 'p', 'T', 'div_U', 'vort', etc.
solution.running_plot.timesteps = 50;                       % Update interval [time steps]

fprintf("Input data processed. \n")
