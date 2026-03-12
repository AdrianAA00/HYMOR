% main.m - Entry point for channel flow stability simulation
%
% Description:
%   Main driver script for running a channel flow stability simulation.
%   This script orchestrates
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
addpath(solver_dir + "chemistry/")
addpath(solver_dir + "utils/Postprocessing/")
addpath(solver_dir + "utils/Time_marching/")
addpath(solver_dir + "utils/Shock_fitting/")
addpath(solver_dir + "utils/Stability_analysis/")
addpath(solver_dir + "utils/Stability_analysis/Eigenvalues/")
addpath(solver_dir + "utils/Stability_analysis/Modal_stability_analysis/")

%% Load input configuration
filename = './input_file.m';
LOAD_INPUT_VARIABLES(filename);

%% Initialize chemistry
chemistry = SET_CHEMISTRY(solution);

%% Initialize mesh, initial condition, and shock position
solution = INITIALIZATION(solution, solution_save, chemistry);

%% Plot initial setup
PLOT_INITIAL_SET_UP(solution);

%% Visualize velocity
figure(1)
hold on
plot_x_coords = solution.mesh.x;
plot_y_coords = solution.mesh.y;
plot_c_data = solution.var.rho_u(2:end-1, 2:end-1) ./ solution.var.rho(2:end-1, 2:end-1);
name_latex = "$u$";
pcolor(plot_x_coords,plot_y_coords,plot_c_data)
cb = colorbar;
cb.Label.String = name_latex;
cb.Label.Interpreter = "latex";
axis equal
shading interp;
colormap(jet);
xlabel("x/h","Interpreter","latex")
ylabel("y/h","Interpreter","latex")
hold off

%% Linearize A
solution.boundary_conditions.periodic = true;
[solution,L] = LINEARIZE_L(solution,chemistry);

%% Compute eigenvalues of stability analysis
n_modes = 4; % search 4 most unstable eigenmodes
[V, D] = eigs(L, n_modes, 0.00 + 0.25i); % Search close to region of reference
for i = 1:n_modes
    fprintf("eigenvalue ")
    fprintf("%i", i)
    fprintf(" = ")
    fprintf("%f + %fi", real(D(i, i)), imag(D(i, i)))
    fprintf("\n")
end

%% Visualize eigenmodes global modal analysis
mode = 1; % Mode to visualize
T_plot = 0; % Time in which disturbance is plotted
freestream_disturbances = false;
PLOT_MODES(freestream_disturbances,L,solution,chemistry,V(:,mode),T_plot); % Plot modes