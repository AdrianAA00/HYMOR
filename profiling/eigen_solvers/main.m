% main.m - Entry point for eigensolver scaling profiling
%
% Description:
%   Test scaling and performance of eigenvalue, transient growth and
%   freestream receptivity solvers across multiple grid sizes.
%   Timings are printed to screen and saved to the "timings/" folder.
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

%% Create timings output directory
timings_dir = fullfile(fileparts(mfilename('fullpath')), 'timings');
if ~exist(timings_dir, 'dir')
    mkdir(timings_dir);
end

%% Load input configuration
filename = './input_file.m';
LOAD_INPUT_VARIABLES(filename);

%% Load chemistry model
if ~solution.restart
    chemistry = SET_CHEMISTRY(solution);
end

%% Run simulation (256 x 64 grid) for base flow
solution.restart        = false;   % Start from previous set up
solution.time_integration.N_iter   = 100;
solution.mesh.Nchi = 2^8;  % Grid points along chi (wall boundary)
solution.mesh.Neta = 2^6;   % Grid points along eta (wall-normal)
solution.freestream.u   = 2000;      % [m/s]
solution.freestream.Re  = 2000;      % Reynolds number (based on tip radius)
solution = INITIALIZATION(solution, solution_save, chemistry);

solution_base = solution;
tic
solution = RUN_SIMULATION(solution, solution_base, chemistry);
toc
solution_save = solution;

% Increase Reynolds for tests
solution.freestream.Re  = 50000;      % Reynolds number (based on tip radius)

%% ========================================================================
%  Grid-scaling timing study
%  ========================================================================

% Define grid sizes to test (Nchi x Neta)
grid_Nchi = [2^8,  2^9,  2^10, 2^11, 2^12];
grid_Neta = [2^6,  2^7,  2^8,  2^9,  2^10];
n_grids   = 3; % Number of grid sizes to test (use first n_grids entries from grid_Nchi and grid_Neta)
n_modes   = 10; % search 10 most unstable eigenmodes
N_l       = 40;

% Pre-allocate timing arrays
t_eigenvalues  = zeros(n_grids, 1);
t_tgd          = zeros(n_grids, 1);
t_freestream   = zeros(n_grids, 1);

% Print header
fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  Eigensolver scaling study — MATLAB\n');
fprintf('==========================================================================\n');
fprintf('%-6s  %-6s  %-12s  %-16s  %-16s  %-16s\n', ...
        'Nchi', 'Neta', 'Nchi*Neta', 'EIGENVALUES [s]', 'TGD [s]', 'FREESTREAM [s]');
fprintf('--------------------------------------------------------------------------\n');

% Write file headers only if files do not exist yet
eig_file = fullfile(timings_dir, 'eigenvalues_timings_MATLAB.txt');
tgd_file = fullfile(timings_dir, 'transient_growth_downstream_timings_MATLAB.txt');
fr_file  = fullfile(timings_dir, 'freestream_receptivity_timings_MATLAB.txt');

if ~exist(eig_file, 'file')
    fid_eig = fopen(eig_file, 'w');
    fprintf(fid_eig, '# Eigensolver scaling: EIGENVALUES\n');
    fprintf(fid_eig, '# %-10s  %-10s  %-14s  %-10s  %-16s\n', 'Nchi', 'Neta', 'Nchi*Neta', 'Language', 'Time [s]');
    fclose(fid_eig);
end
if ~exist(tgd_file, 'file')
    fid_tgd = fopen(tgd_file, 'w');
    fprintf(fid_tgd, '# Eigensolver scaling: TRANSIENT_GROWTH_DOWNSTREAM\n');
    fprintf(fid_tgd, '# %-10s  %-10s  %-14s  %-10s  %-16s\n', 'Nchi', 'Neta', 'Nchi*Neta', 'Language', 'Time [s]');
    fclose(fid_tgd);
end
if ~exist(fr_file, 'file')
    fid_fr = fopen(fr_file, 'w');
    fprintf(fid_fr, '# Eigensolver scaling: FREESTREAM_RECEPTIVITY\n');
    fprintf(fid_fr, '# %-10s  %-10s  %-14s  %-10s  %-16s\n', 'Nchi', 'Neta', 'Nchi*Neta', 'Language', 'Time [s]');
    fclose(fid_fr);
end

for ig = 1:n_grids

    Nchi_i = grid_Nchi(ig);
    Neta_i = grid_Neta(ig);

    %% Set up grid and linearize
    solution.mesh.Nchi = Nchi_i;
    solution.mesh.Neta = Neta_i;
    solution.restart   = true;
    solution = INITIALIZATION(solution, solution_save, chemistry);
    [solution, L] = LINEARIZE_L(solution, chemistry);

    %% Time EIGENVALUES (use tic handle to avoid confusion with inner tic/toc)
    tStart_eig = tic;
    [~, ~] = EIGENVALUES(L, solution, n_modes);
    t_eigenvalues(ig) = toc(tStart_eig);

    %% Time TRANSIENT_GROWTH_DOWNSTREAM
    T_TGD = 1;
    V_TGD = zeros(4 * Nchi_i * Neta_i, n_modes, size(T_TGD,1));
    D_TGD = zeros(n_modes, n_modes, length(T_TGD));

    tStart_tgd = tic;
    [~, ~, ~] = TRANSIENT_GROWTH_DOWNSTREAM(L, solution, n_modes, T_TGD);
    t_tgd(ig) = toc(tStart_tgd);

    %% Linearize extended system
    w_infty = 2*pi*linspace(0, N_l, N_l+1)';
    [solution, L_] = LINEARIZE_L_(L, solution, chemistry, w_infty);

    %% Time FREESTREAM_RECEPTIVITY
    T_TGF = 1;
    V_TGF = zeros(4 * Nchi_i * Neta_i + 4 * Nchi_i * (N_l + 1), n_modes, length(T_TGF));
    D_TGF = zeros(n_modes, n_modes, length(T_TGF));

    tStart_fr = tic;
    [~, ~, ~] = FREESTREAM_RECEPTIVITY(L_, solution, n_modes, T_TGF, w_infty);
    t_freestream(ig) = toc(tStart_fr);

    %% Print row
    fprintf('%-6d  %-6d  %-12d  %-16.4f  %-16.4f  %-16.4f\n', ...
            Nchi_i, Neta_i, Nchi_i*Neta_i, ...
            t_eigenvalues(ig), t_tgd(ig), t_freestream(ig));

    %% Append this iteration's results to timing files
    fid_eig = fopen(eig_file, 'a');
    fprintf(fid_eig, '  %-10d  %-10d  %-14d  %-10s  %-16.6f\n', ...
            Nchi_i, Neta_i, Nchi_i*Neta_i, 'MATLAB', t_eigenvalues(ig));
    fclose(fid_eig);

    fid_tgd = fopen(tgd_file, 'a');
    fprintf(fid_tgd, '  %-10d  %-10d  %-14d  %-10s  %-16.6f\n', ...
            Nchi_i, Neta_i, Nchi_i*Neta_i, 'MATLAB', t_tgd(ig));
    fclose(fid_tgd);

    fid_fr = fopen(fr_file, 'a');
    fprintf(fid_fr, '  %-10d  %-10d  %-14d  %-10s  %-16.6f\n', ...
            Nchi_i, Neta_i, Nchi_i*Neta_i, 'MATLAB', t_freestream(ig));
    fclose(fid_fr);

    %% Free memory
    clear L L_

end

fprintf('==========================================================================\n\n');
fprintf('Timing files saved in: %s\n', timings_dir);