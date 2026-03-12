% main.m - Entry point for 3D axisymmetric blunted cone simulation
%
% Description:
%   Main driver script for running a hypersonic flow stability simulation
%   over a 3D axisymmetric blunted cone geometry. This script orchestrates
%   the full simulation pipeline:
%     1. Adds required solver utility paths
%     2. Loads user-defined input parameters from input_file.m
%     3. Initializes the chemistry model (if not restarting)
%     4. Sets up the computational mesh and initial conditions
%     5. Runs the time-marching simulation
%     6. Generates post-processing plots
%
% Inputs:
%   None (all parameters are loaded from input_file.m)
%
% Outputs:
%   solution     - struct : Final solution state containing flow fields,
%                           mesh data, and solver metadata
%   chemistry    - struct : Chemistry model data (species, reactions, etc.)
%
% Dependencies:
%   - input_file.m (input configuration file in same directory)
%   - Solver utility modules:
%       utils/Initialization/    - LOAD_INPUT_VARIABLES, INITIALIZATION
%       utils/Mesh/              - Mesh generation and adaptation
%       utils/Operators/         - Discrete differential operators
%       utils/Energy_budgets/    - Energy budget analysis tools
%       utils/Postprocessing/    - PLOT_INITIAL_SET_UP, PLOTS
%       utils/Shock_fitting/     - Shock fitting routines
%       utils/Time_marching/     - RUN_SIMULATION, time integrators
%       chemistry/               - SET_CHEMISTRY, thermochemical models
%
% Usage:
%   Set solver_dir below to the absolute path of the solver root directory,
%   then run this script from the test case directory:
%       >> run('main.m')
%
% Notes:
%   - The solver_dir path must be updated to match the local installation.
%
% Part of: Hypersonics Stability MATLAB Solver - Test Cases Module

%% Parent directory
clear all
clc
solver_dir = '../../'; % Update this path to the absolute path of the solver root directory on your system
solution = struct();
solution.solver_dir = solver_dir;

%% Include utilities
addpath(solver_dir + "utils/Initialization/")
addpath(solver_dir + "utils/Mesh/")
addpath(solver_dir + "utils/Operators/")
addpath(solver_dir + "utils/Energy_budgets/")
addpath(solver_dir + "utils/Postprocessing/")
addpath(solver_dir + "utils/Shock_fitting/")
addpath(solver_dir + "utils/Time_marching/")
addpath(solver_dir + "chemistry/")
addpath(solver_dir + "utils/Stability_analysis/")
addpath(solver_dir + "utils/Stability_analysis/Eigenvalues/")
addpath(solver_dir + "utils/Stability_analysis/Modal_stability_analysis/")
addpath(solver_dir + "utils/Stability_analysis/Transient_growth_downstream/")
addpath(solver_dir + "utils/Stability_analysis/Freestream_receptivity/")

%% Load input configuration
filename = './input_file.m';
LOAD_INPUT_VARIABLES(filename);

%% Load chemistry model
if ~solution.restart
    chemistry = SET_CHEMISTRY(solution);
end

%% Plot initial setup
solution = INITIALIZATION(solution, solution_save, chemistry);
PLOT_INITIAL_SET_UP(solution);

%% Run simulation (coarse grid lower speed for convergence)
solution.restart        = false;      % Start from scratch
solution.time_integration.N_iter   = 1500;
solution.freestream.u   = 2000;       % [m/s]
solution.freestream.Re  = 2000;      % Reynolds number (based on tip radius)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Run simulation (coarse grid increase speed)
% Need to increase gradually freestream speed, because if not strong shock can be
% generated after the main shock and can break numerics
solution.restart        = true;      % Restart from previous case
solution.shock.interpolate              = "1st"; % Set first order interpolation for increased robustness in strong shock transient
solution.time_integration.N_iter   = 3500;
solution.freestream.u   = 5000;     % [m/s]
solution.freestream.Re  = 2000;     % Reynolds number (based on tip radius)
solution = INITIALIZATION(solution, solution_save, chemistry); % Restart from previous case

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Run simulation (remesh to finer grid & higher Reynolds)
solution.restart        = true;      % Restart from previous case
solution.shock.interpolate              = "2nd";
solution.time_integration.N_iter   = 1000; % Increase more iterations in case want better converged case
solution.freestream.u   = 5000;      % [m/s]
solution.freestream.Re  = 20000;     % Reynolds number (based on tip radius)
solution.mesh.Nchi = 400;  % Points along chi (wall boundary)
solution.mesh.Neta = 80;   % Points along eta (wall-normal direction)
solution.curvilinear_mapping.refinement_stagnation.state        = true;
solution.curvilinear_mapping.refinement_stagnation.BL_thickness = 0.2; % Thickness region that is refined
solution.curvilinear_mapping.refinement_stagnation.intensity    = 0.95; % 0 - not refined, 1 - infinitely refined
solution.curvilinear_mapping.refinement_wall.state              = true;
solution.curvilinear_mapping.refinement_wall.BL_thickness       = 0.1; % Thickness region that is refined
solution.curvilinear_mapping.refinement_wall.intensity          = 0.995; % 0 - not refined, 1 - infinitely refined
solution = INITIALIZATION(solution, solution_save, chemistry); % Restart from previous case

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Visualize density and vorticity
figure(1)
hold on
plot_x_coords = solution.mesh.x;
plot_y_coords = solution.mesh.y;
plot_c_data = solution.var.rho(2:end-1, 2:end-1) ./ solution.freestream.rho_0;
name_latex = "$\rho/\rho_\infty$";
pcolor(plot_x_coords,plot_y_coords,plot_c_data)
cb = colorbar;
cb.Label.String = name_latex;
cb.Label.Interpreter = "latex";
axis equal
shading interp;
colormap(jet);
xlabel("x/R","Interpreter","latex")
ylabel("y/R","Interpreter","latex")
hold off

figure(2)
hold on
plot_x_coords = solution.mesh.x(2:end-1, 2:end-1);
plot_y_coords = solution.mesh.y(2:end-1, 2:end-1);
u_vel_int = solution.var.rho_u(2:end-1, 2:end-1) ./ solution.var.rho(2:end-1, 2:end-1);
v_vel_int = solution.var.rho_v(2:end-1, 2:end-1) ./ solution.var.rho(2:end-1, 2:end-1);
[~, d_u_dy] = DERIVATIVE(u_vel_int, solution);
[d_v_dx, ~] = DERIVATIVE(v_vel_int, solution);
plot_c_data = d_v_dx - d_u_dy;
plot_c_data = plot_c_data / solution.freestream.U * solution.curvilinear_mapping.R;
name_latex = "$\omega R/U_\infty$";
pcolor(plot_x_coords, plot_y_coords, plot_c_data)
cb = colorbar;
cb.Label.String = name_latex;
cb.Label.Interpreter = "latex";
axis equal
xlabel("x/R","Interpreter","latex")
ylabel("y/R","Interpreter","latex")
shading interp;
colormap(jet);
hold off

%% Linearize A
[solution,L] = LINEARIZE_L(solution,chemistry);

%% Compute eigenvalues of stability analysis
n_modes = 10; % search 10 most unstable eigenmodes
[V,D] = EIGENVALUES(L,solution,n_modes); % V stores eigenmodes, and D the eigenvalues

%% Visualize eigenmodes global modal analysis
mode = 3; % Mode to visualize
T_plot = 0; % Time in which disturbance is plotted
freestream_disturbances = false;
PLOT_MODES(freestream_disturbances,L,solution,chemistry,V(:,mode),T_plot); % Plot modes

%% Compute transient growth downstream flow
T_TGD = [0.5;1;1.5]; % Time to find optimal transient growth
n_modes = 5; % Number of modes to find in each iteration
V_TGD = zeros(4 * solution.mesh.Nchi * solution.mesh.Neta, n_modes, size(T_TGD,1)); % Store modes
D_TGD = zeros(n_modes, n_modes, length(T_TGD)); % Store transient growth

for i = 1:length(T_TGD)
    [V_TGD(:,:,i),D_TGD(:,:,i),T_opt_TGD(i,1)] = TRANSIENT_GROWTH_DOWNSTREAM(L,solution,n_modes,T_TGD(i,1));
end

%% Visualize eigenmodes  global non-modal analysis
time_optimization_index = 2;
mode = 1; % Mode to visualize
T_plot = 1.50; % Time in which disturbance is plotted
freestream_disturbances = false;
PLOT_MODES(freestream_disturbances,L,solution,chemistry,V_TGD(:,mode,time_optimization_index),T_plot); % Plot modes

%% Linearize A_ (extended system with freestream perturbations)
N_l = 40;
w_infty = 2*pi*linspace(0,N_l,N_l+1)'; % Streamwise frequencies allowed of freestream disturbances

[solution,L_] = LINEARIZE_L_(L,solution,chemistry,w_infty);

%% Compute freestream receptivity optimal modes
T_TGF = [5]; % Time to find optimal transient growth
n_modes = 5;
V_TGF = zeros(4 * solution.mesh.Nchi * solution.mesh.Neta + 4 * solution.mesh.Nchi * (N_l + 1), n_modes, length(T_TGF));
D_TGF = zeros(n_modes, n_modes, length(T_TGF));

for i = 1:length(T_TGF)
    [V_TGF(:,:,i),D_TGF(:,:,i),T_opt_TGF(i,1)] = FREESTREAM_RECEPTIVITY(L_,solution,n_modes,T_TGF(i,1),w_infty);
end

%% Visualize freestream receptivity modes
time_optimization_index = 1;
mode = 2; % Mode to visualize
T_plot = 5; % Time in which disturbance is plotted
freestream_disturbances = true;
PLOT_MODES(freestream_disturbances,L_,solution,chemistry,V_TGF(:,mode,time_optimization_index),T_plot,w_infty); % Plot modes
