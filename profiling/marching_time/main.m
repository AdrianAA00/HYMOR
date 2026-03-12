% main.m - Entry point for 3D axisymmetric blunted cone simulation
%
% Description:
%   Script to measure scaling of the time-marching solver without chemistry and no shock-fitting
%
% Usage:
%   Set solver_dir below to the absolute path of the solver root directory,
%   then run this script from the test case directory:
%       >> run('main.m')
%
% Notes:
%   - The solver_dir path must be updated to match the local installation.
%
% Part of: Hypersonics Stability MATLAB Solver - Profiling Module

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
addpath(solver_dir + "utils/Shock_fitting/")
addpath(solver_dir + "utils/Postprocessing/")
addpath(solver_dir + "utils/Time_marching/")

%% Load input configuration
filename = './input_file.m';
LOAD_INPUT_VARIABLES(filename);

%% Initialize chemistry
chemistry = SET_CHEMISTRY(solution);

%% Open output file for timing results
results_file = fopen('./timings/timing_results_matlab.txt', 'w');
fprintf(results_file, '# Timing results - MATLAB\n');
fprintf(results_file, '# Grid (Nchi x Neta)    Elapsed time (s)\n');
fprintf(results_file, '%-25s %s\n', 'Grid', 'Elapsed_time_s');

%% Time marching 100 x 100 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 100;  % Grid points along chi (streamwise)
solution.mesh.Neta = 100;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

%% Time marching 200 x 200 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 200;  % Grid points along chi (streamwise)
solution.mesh.Neta = 200;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

%% Time marching 400 x 400 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 400;  % Grid points along chi (streamwise)
solution.mesh.Neta = 400;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

%% Time marching 800 x 800 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 800;  % Grid points along chi (streamwise)
solution.mesh.Neta = 800;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

%% Time marching 1600 x 1600 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 1600;  % Grid points along chi (streamwise)
solution.mesh.Neta = 1600;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

%% Time marching 3200 x 3200 grid, 1000 iterations
solution.restart                    = false;
solution.time_integration.N_iter    = 1000;          % Number of time steps
solution.mesh.Nchi = 1600;  % Grid points along chi (streamwise)
solution.mesh.Neta = 1600;  % Grid points along eta (wall-normal)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;

tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
elapsed = toc;
fprintf('Grid: %d x %d - Elapsed time: %.6f s\n', solution.mesh.Nchi, solution.mesh.Neta, elapsed);
fprintf(results_file, '%-25s %.6f\n', sprintf('%dx%d', solution.mesh.Nchi, solution.mesh.Neta), elapsed);
clear solution_base

fclose(results_file);
