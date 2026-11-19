% main.m - Non-axisymmetric stability of Hagen-Poiseuille pipe flow
%
% Description:
%   Modal and non-modal stability of laminar pipe flow to disturbances
%   with azimuthal wavenumber m, computed on the axisymmetric (r, x) grid
%   with stability_analysis.azimuthal_wavenumber. The script
%     1. Loads the inputs and initializes the Hagen-Poiseuille base flow
%     2. Sweeps m = 0..3: linearizes for each m and computes the least
%        stable eigenvalues of axial wavenumber alpha = 2*pi/Lx
%     3. Computes the optimal transient growth of streamwise-invariant
%        (alpha = 0) disturbances with m = 1 (lift-up)
%   and compares every result with the spectral reference solution of
%   PIPE_FLOW_REFERENCE (incompressible, Chebyshev collocation).
%
% Inputs:
%   None (all parameters are loaded from input_file.m)
%
% Outputs:
%   modal     - struct array : HYMOR and reference eigenvalues per m
%   nonmodal  - struct       : HYMOR and reference optimal gains G(T)
%
% Dependencies:
%   - input_file.m, PIPE_FLOW_REFERENCE.m (this directory)
%   - Solver utility modules (utils/...)
%
% Part of: Hypersonics Stability MATLAB Solver - Tutorials (pipe flow)

clear all
clc

solver_dir = '../../'; % Update this path to the absolute path of the solver root directory on your system
solution = struct();
solution.solver_dir = solver_dir;

%% Add paths
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
addpath(solver_dir + "utils/Stability_analysis/Transient_growth_downstream/")

%% Load input data and initialize the base flow
filename = './input_file.m';
LOAD_INPUT_VARIABLES(filename);
chemistry = SET_CHEMISTRY(solution);
solution = INITIALIZATION(solution, solution_save, chemistry);

Re = solution.freestream.Re;
alpha = 2*pi / solution.curvilinear_mapping.Lx;

%% 1. Modal analysis: parametric sweep over the azimuthal wavenumber
m_sweep = 0:3;   % Azimuthal wavenumbers
n_modes = 3;     % Least stable modes compared per m
modal = struct('m', num2cell(m_sweep), 'lambda', [], 'lambda_ref', []);

for k = 1:numel(m_sweep)
    solution.stability_analysis.azimuthal_wavenumber = m_sweep(k);
    [~, L] = LINEARIZE_L(solution, chemistry);

    % exp(i alpha x) disturbances: HYMOR's real operator holds both signs of
    % alpha, so the reference eigenvalues are searched as complex conjugates
    reference = PIPE_FLOW_REFERENCE(Re, alpha, m_sweep(k), 60);
    modal(k).lambda_ref = conj(reference.lambda(1:n_modes));
    modal(k).lambda = arrayfun(@(lr) eigs(L, 1, lr), modal(k).lambda_ref);
end

fprintf("\nLeast stable eigenvalues, alpha = %.2f, Re = %g (HYMOR | reference | relative error)\n", alpha, Re);
for k = 1:numel(m_sweep)
    for j = 1:n_modes
        fprintf("  m = %d:  %+.5f %+.5fi | %+.5f %+.5fi | %.1e\n", modal(k).m, ...
            real(modal(k).lambda(j)), imag(modal(k).lambda(j)), ...
            real(modal(k).lambda_ref(j)), imag(modal(k).lambda_ref(j)), ...
            abs(modal(k).lambda(j) - modal(k).lambda_ref(j)) / abs(modal(k).lambda_ref(j)));
    end
end

figure()
hold on
colors = lines(numel(m_sweep));
for k = 1:numel(m_sweep)
    plot(imag(modal(k).lambda_ref), real(modal(k).lambda_ref), 'o', 'Color', colors(k, :), ...
        'MarkerSize', 9, 'DisplayName', sprintf('reference, $m = %d$', modal(k).m))
    plot(imag(modal(k).lambda), real(modal(k).lambda), 'x', 'Color', colors(k, :), ...
        'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', sprintf('HYMOR, $m = %d$', modal(k).m))
end
xlabel("$\omega = \Im(\lambda)$", "Interpreter", "latex")
ylabel("$\sigma = \Re(\lambda)$", "Interpreter", "latex")
legend("Interpreter", "latex", "Location", "best")
box on
hold off

%% 2. Non-modal analysis: lift-up of streamwise-invariant disturbances (alpha = 0, m = 1)
% Three cells in x with a short period resolve only alpha = 0 (a longer
% period would also carry unresolved alpha ~= 0 content). The time step of
% the exponential time-stepper is acoustic, dt ~ Mach: Mach = 0.1 halves the
% cost and changes G by less than 0.2%.
solution_tg = solution;
solution_tg.mesh.Nchi = 3;
solution_tg.mesh.Neta = 24;
solution_tg.curvilinear_mapping.Lx = 0.1;
solution_tg.freestream.Mach = 0.1;
solution_tg = INITIALIZATION(solution_tg, [], chemistry);
solution_tg.stability_analysis.azimuthal_wavenumber = 1;
[solution_tg, L_tg] = LINEARIZE_L(solution_tg, chemistry);

nonmodal.T = [10; 50];          % Time horizons
nonmodal.G = zeros(size(nonmodal.T));
for i = 1:numel(nonmodal.T)
    [~, D_tg, nonmodal.T(i)] = TRANSIENT_GROWTH_DOWNSTREAM(L_tg, solution_tg, 1, nonmodal.T(i));
    nonmodal.G(i) = real(D_tg(1, 1));
end
reference_tg = PIPE_FLOW_REFERENCE(Re, 0, 1, 60);
nonmodal.G_ref = reference_tg.gain(nonmodal.T);

fprintf("\nOptimal energy gain, alpha = 0, m = 1, Re = %g (HYMOR | reference)\n", Re);
for i = 1:numel(nonmodal.T)
    fprintf("  T = %5.1f:  G = %8.3f | %8.3f\n", nonmodal.T(i), nonmodal.G(i), nonmodal.G_ref(i));
end

figure()
T_curve = linspace(0, 100, 201);
plot(T_curve, reference_tg.gain(T_curve), 'k-', 'DisplayName', 'reference')
hold on
plot(nonmodal.T, nonmodal.G, 'rx', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'HYMOR')
xlabel("$T U_c / R$", "Interpreter", "latex")
ylabel("$G(T)$", "Interpreter", "latex")
legend("Interpreter", "latex", "Location", "best")
box on
hold off
