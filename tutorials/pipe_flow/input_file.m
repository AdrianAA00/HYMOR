% input_file.m - Input configuration for Hagen-Poiseuille pipe flow
%
% Laminar flow in a pipe of radius R = 1 (axis at y = 0, wall at y = R),
% periodic in the axial direction x over one wavelength 2*pi/alpha. The
% base flow is the Hagen-Poiseuille profile u = 1 - (r/R)^2 (centreline
% velocity scale, Re = U_c R / nu) at low Mach number, linearized for
% non-axisymmetric disturbances with azimuthal wavenumber m.
%
% Part of: Hypersonics Stability MATLAB Solver - Tutorials (pipe flow)

%% Restart options
solution.restart               = false;  % Resume from a previous run
solution.remesh                = false;  % Regenerate the mesh on restart

if ~exist('solution_save', 'var') || isempty(solution_save)
    solution_save = [];
end

%% Governing equations
solution.PDE_dimension = "3D-axisymmetric";  % Options: "2D", "3D-axisymmetric"

%% Chemistry
solution.chemistry.is_chemistry_enabled = false;  % Calorically perfect gas

%% Mesh
solution.mesh.Nchi = 48;  % Cells along chi (axial, multiple of 3 for periodic BCs)
solution.mesh.Neta = 64;  % Cells along eta (radial, axis to wall)

%% Curvilinear mapping
solution.curvilinear_mapping.boundary_type        = "pipe";  % Geometry type
solution.curvilinear_mapping.Lx                   = 2*pi;    % Axial period (alpha = 2*pi/Lx = 1)
solution.curvilinear_mapping.R                    = 1;       % Pipe radius
solution.curvilinear_mapping.eta_refinement_power = 1.3;     % Radial refinement towards the wall (1 = uniform)

%% Boundary conditions
solution.boundary_conditions.boundary_eta0.name = 'symmetry';           % Pipe axis
solution.boundary_conditions.boundary_eta1.name = 'no_slip_adiabatic';  % Pipe wall
solution.boundary_conditions.boundary_chi0.name = 'periodic';           % Axial periodicity
solution.boundary_conditions.boundary_chi1.name = 'periodic';

%% Freestream conditions (centreline values)
solution.freestream.u     = 1;       % Centreline velocity
solution.freestream.v     = 0;       % Radial velocity
solution.freestream.gamma = 1.4;     % Specific heat ratio
solution.freestream.Mach  = 0.05;    % Centreline Mach number (nearly incompressible)
solution.freestream.Pr    = 0.73;    % Prandtl number
solution.freestream.Re    = 1000;    % Reynolds number U_c R / nu

%% Time integration (sets the step of the exponential time-stepper)
solution.time_integration.N_iter          = 1;              % No base-flow time marching
solution.time_integration.time_integrator = "Explicit_RK4"; % Scheme: "Explicit_RK4", "Implicit_Euler"
solution.time_integration.CFL             = 1;              % CFL number of the linear time-stepper
solution.time_integration.dt              = 1e-6;           % Initial time step
solution.time_integration.max_dt          = 0.3;            % Maximum allowable time step

%% Shock fitting
solution.shock.enabled = false;  % No shock

%% Stability analysis
solution.stability_analysis.perturbation_magnitude = 1e-6;      % Linearization perturbation size (rho*E ~ 1/Mach^2: 1e-8 loses digits)
solution.stability_analysis.eigenvalue_solver      = "CPU_LU";  % Eigensolver: "CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"
solution.stability_analysis.device                 = "auto";    % Time-stepper device: "auto" (GPU for large operators), "cpu", "gpu"
solution.stability_analysis.azimuthal_wavenumber   = 1;         % Azimuthal wavenumber m of the disturbances

%% Running plots
solution.running_plot.enabled = false;  % Show live plots during simulation

fprintf("Input data processed. \n")
