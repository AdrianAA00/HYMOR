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
%   solution_old - struct : Solution from the previous time step
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

%% 
PLOT_INITIAL_SET_UP(solution);

%% Run simulation (coarse grid for fast convergence)
solution.restart        = false;     % Start from scratch
solution.time_integration.N_iter   = 2000;
solution.freestream.Re  = 10000;     % Reynolds number 
solution.shock.spline_param        = 1 - 1e-6 * 100^3 / solution.mesh.Nchi^3; % Smoothing parameter for csaps
solution.mesh.Nchi = 200;  % Points along chi (wall boundary)
solution.mesh.Neta = 40;   % Points along eta (wall-normal direction)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Run simulation (fine grid)
solution.restart        = true;     % Start from scratch
solution.time_integration.N_iter   = 1000;
solution.freestream.Re  = 100000;     % Reynolds number
solution.shock.spline_param        = 1 - 1e-8 * 100^3 / solution.mesh.Nchi^3; % Smoothing parameter for csaps
solution.mesh.Nchi = 800;  % Points along chi (wall boundary)
solution.mesh.Neta = 100;   % Points along eta (wall-normal direction)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Run simulation (finest grid) (if this refinement it takes 2 - 3 hours to run)
solution.restart        = true;     % Start from scratch
solution.time_integration.N_iter   = 1000;
solution.freestream.Re  = 100000;     % Reynolds number
solution.mesh.Nchi = 1600;  % Points along chi (wall boundary)
solution.mesh.Neta = 200;   % Points along eta (wall-normal direction)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

%% Visualize density and vorticity
label_size = 20;
legend_size = 20;
tick_size = 15;

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
cb.Label.FontSize = legend_size;
cb.TickLabelInterpreter = 'latex';
axis equal
shading interp;
colormap(jet);
xlabel("x/R","Interpreter","latex", "FontSize", label_size)
ylabel("y/R","Interpreter","latex", "FontSize", label_size)
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', tick_size)
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
cb.Label.FontSize = legend_size;
cb.TickLabelInterpreter = 'latex';
axis equal
xlabel("x/R","Interpreter","latex", "FontSize", label_size)
ylabel("y/R","Interpreter","latex", "FontSize", label_size)
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', tick_size)
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
solution.running_plot.scaling_range = 1/2; % Saturate to see field
PLOT_MODES(freestream_disturbances,L,solution,chemistry,V(:,mode),T_plot); % Plot modes

%% Compute eigenvalues of transient growth downstream flow
T_TGD = [3]; % Time to find optimal transient growth
n_modes = 5; % Number of modes to find in each iteration
V_TGD = zeros(4 * solution.mesh.Nchi * solution.mesh.Neta, n_modes, size(T_TGD,1)); % Store modes
D_TGD = zeros(n_modes, n_modes, length(T_TGD)); % Store transient growth

for i = 1:length(T_TGD)
    [V_TGD(:,:,i),D_TGD(:,:,i),T_opt_TGD(i,1)] = TRANSIENT_GROWTH_DOWNSTREAM(L,solution,n_modes,T_TGD(i,1));
end

%% Visualize eigenmodes  global non-modal analysis
time_optimization_index = 1;
mode = 2; % Mode to visualize (the first mode is a spurious numerical mode, reflection at the boundary)
T_plot = 0.0; % Time in which disturbance is plotted
freestream_disturbances = false;
solution.running_plot.scaling_range = 1/2; % Saturate to see field
PLOT_MODES(freestream_disturbances,L,solution,chemistry,V_TGD(:,mode,time_optimization_index),T_plot); % Plot modes

%% Visualize selected non-modal
T_f = 5; % Final time of integration to see energy growth
mode = 2; % Mode to evolve in time
freestream_disturbances = false;
get_amplification_only = true;
max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances,V_TGD(:,mode,time_optimization_index),D_TGD(mode,mode,time_optimization_index),T_opt_TGD(time_optimization_index,1),L,solution,chemistry,T_f,get_amplification_only);

disp("max_ref(E) = "   + max_gain.non_temporal.all)
disp("max_ref(E_p) = " + max_gain.non_temporal.acoustic)
disp("max_ref(E_S) = " + max_gain.non_temporal.entropic)
disp("max_ref(E_k) = " + max_gain.non_temporal.kinetic)

%% Linearize A_ (extended system with freestream perturbations)
N_l = 40;
w_infty = 2*pi*linspace(0,N_l,N_l+1)'; % Streamwise frequencies allowed of freestream disturbances
[solution,L_] = LINEARIZE_L_(L,solution,chemistry,w_infty);

%% Compute optimal freestream receptivity modes
T_TGF = [5]; % Time to find optimal transient growth
n_modes = 5;
V_TGF = zeros(4 * solution.mesh.Nchi * solution.mesh.Neta + 4 * solution.mesh.Nchi * (N_l + 1), n_modes, length(T_TGF));
D_TGF = zeros(n_modes, n_modes, length(T_TGF));

for i = 1:length(T_TGF)
    [V_TGF(:,:,i),D_TGF(:,:,i),T_opt_TGF(i,1)] = FREESTREAM_RECEPTIVITY(L_,solution,n_modes,T_TGF(i,1),w_infty);
end

%% Visualize freestream receptivity modes
time_optimization_index = 1;
mode = 3; % Mode to visualize (first two spurious)
T_plot = 10; % Time in which disturbance is plotted
freestream_disturbances = true;
solution.running_plot.scaling_range = 1/10; % Saturate to see field
PLOT_MODES(freestream_disturbances,L_,solution,chemistry,V_TGF(:,mode,time_optimization_index),T_plot,w_infty); % Plot modes

%% Visualize selected freestrean gains
T_f = 10; % Final time of integration to see energy growth
mode = 3; % Mode to evolve in time
freestream_disturbances = true;
get_amplification_only = true;
max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances,V_TGF(:,mode,time_optimization_index),D_TGF(mode,mode,time_optimization_index),T_opt_TGF(time_optimization_index,1),L_,solution,chemistry,T_f,get_amplification_only,w_infty);

disp("max_ref(E) = "   + max_gain.non_temporal.all)
disp("max_ref(E_p) = " + max_gain.non_temporal.acoustic)
disp("max_ref(E_S) = " + max_gain.non_temporal.entropic)
disp("max_ref(E_k) = " + max_gain.non_temporal.kinetic)

%% Save workspace
save('myWorkspace.mat')
