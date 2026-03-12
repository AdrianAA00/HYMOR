% input_file.m - Input configuration for channel flow
%
% Description:
%   Defines all simulation parameters for a subsonic channel flow case.
%   This configuration file is loaded by main.m via LOAD_INPUT_VARIABLES
%   and populates the 'solution' struct.
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

solution.restart               = false;  % Resume from a previous run
solution.remesh                = false;  % Regenerate the mesh on restart

if ~exist('solution_save', 'var') || isempty(solution_save)
    solution_save = [];
end

%% ========================================================================
%  Simulation Parameters
%  ========================================================================

%% PDE dimension
solution.PDE_dimension = "2D";  % Options: "2D", "3D-axisymmetric"

%% Chemistry model
solution.chemistry.is_chemistry_enabled = false;  % Enable thermochemistry (false = calorically perfect gas)

%% Mesh parameters
solution.mesh.Nchi = 180;  % Grid points along chi (streamwise)
solution.mesh.Neta = 60;   % Grid points along eta (wall-normal)

%% Geometry: channel flow
solution.curvilinear_mapping.boundary_type = "channel";       % Geometry type
solution.curvilinear_mapping.channel_angle = 0 * pi / 8;     % Channel inclination angle [rad]
solution.curvilinear_mapping.eta_refinement_power = 1.6;      % Wall-normal refinement (1 = uniform, >1 = refined near wall)
solution.curvilinear_mapping.Lx       = 2*pi;                % Streamwise domain length
solution.curvilinear_mapping.Ly       = 2;                   % Wall-normal domain height

%% Boundary conditions
%         _____eta1_______
%        |                |
%   chi0 |                | chi1
%        |                |
%        ------eta0--------
%
% Available types:
%   'inflow_subsonic'  'inflow_supersonic'  'periodic'  'shock'
%   'outflow_supersonic'  'outflow_subsonic' 'outflow_NRCBC'
%   'no_slip_adiabatic'  'no_slip_isothermal'  'symmetry'

solution.boundary_conditions.boundary_eta0.name        = 'no_slip_adiabatic';  % Bottom wall
solution.boundary_conditions.boundary_eta1.name        = 'no_slip_adiabatic';  % Top wall
solution.boundary_conditions.boundary_chi0.name        = 'periodic';           % Inflow
solution.boundary_conditions.boundary_chi1.name        = 'periodic';           % Outflow

%% Freestream conditions (non-chemistry mode)
solution.freestream.u       = 1;       % Freestream velocity
solution.freestream.v       = 0;       % Transverse velocity
solution.freestream.gamma   = 1.4;     % Specific heat ratio
solution.freestream.Mach    = 0.05;    % Freestream Mach number
solution.freestream.Pr      = 0.73;    % Prandtl number
solution.freestream.Re      = 7500;    % Reynolds number

%% Time integration
solution.time_integration.N_iter            = 1000;          % Number of time steps
solution.time_integration.time_integrator   = "Explicit_RK4"; % Scheme: "Explicit_RK4", "Implicit_Euler"
solution.time_integration.CFL               = 2;             % CFL number for adaptive time stepping
solution.time_integration.dt                = 0.000001;      % Initial time step [s]
solution.time_integration.max_dt            = 0.3;           % Maximum allowable time step [s]

%% Shock fitting
solution.shock.enabled  = false;  % Enable shock fitting

%% Stability analysis and transient growth
solution.stability_analysis.perturbation_magnitude  = 1e-8;                    % Linearization perturbation size
solution.stability_analysis.eigenvalue_solver       = "GPU_TIMESTEPPER_ARNOLDI"; % Eigensolver: "CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"

%% Live plots during time-marching
solution.running_plot.enabled   = false;                    % Show live plots during simulation

fprintf("Input data processed. \n")
