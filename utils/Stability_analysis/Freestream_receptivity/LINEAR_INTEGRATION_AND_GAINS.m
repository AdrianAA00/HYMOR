function max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances, v, lambda, T_opt, A_, s, chemistry, T_f, get_amplification_only, w_infty)
% LINEAR_INTEGRATION_AND_GAINS  Integrate the linearized system, plot energy gains, and compute budgets.
%
%   max_gain = LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances, v, lambda, T_opt, A_, s, chemistry, T_f, w_infty, get_amplification_only)
%
%   Propagates the optimal initial perturbation v forward in time under the
%   linearized dynamics defined by A_ using a GPU-accelerated Taylor-series
%   matrix exponential. Computes the decomposed energy evolution (total,
%   pressure, kinetic, entropy) and generates publication-quality figures
%   of energy gains and budgets.
%
%   Inputs:
%       freestream_disturbances - If true, use the extended system (downstream + freestream)
%       v                      - Optimal initial perturbation eigenvector
%       lambda                 - Optimal gain (eigenvalue) for overlay on plots
%       T_opt                  - Optimal time horizon corresponding to lambda
%       A_                     - System matrix (or extended system matrix)
%       s               - Solution structure with flow fields and mesh data
%       chemistry              - Chemistry model structure
%       T_f                    - Final integration time for plotting
%       w_infty                - Vector of freestream disturbance frequencies
%       get_amplification_only - If true, skip budget computations and sampling
%
%   Outputs:
%       max_gain - Structure with fields:
%                    .temporal.all, .acoustic, .kinetic, .entropic
%                    .non_temporal.all, .acoustic, .kinetic, .entropic
%                    .budgets (only when ~get_amplification_only)
%
%   Notes:
%       - Uses 5th-order Taylor expansion for the matrix exponential.
%       - All heavy computations are performed on the GPU.
%       - Generates multiple publication-quality figures for energy gains,
%         kinetic energy budgets, and entropic energy budgets.
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Freestream Module

    %% Time-stepping parameters
    orderTaylor = 5; % Order of Taylor expansion for matrix exponential approximation

    s = CFL_TIMESTEP(s);
    t_temp = s.time_integration.dt;
    n_t = round(T_f / t_temp);
    dt = T_f / n_t;
    N_skip = s.mesh.Neta / 2;
    N_samples = max(floor(n_t / N_skip), 500);

    fprintf("\n")
    fprintf("Number timesteps = %f ", n_t)
    fprintf("\n")

    %% Construct energy-norm operators
    if freestream_disturbances
        R_ = CONSTRUCT_R_(s, w_infty);
        norms = [1, 1, 1, 1];
        M_ = CONSTRUCT_M_(s, norms, w_infty);

        T = 1;
        scaling_non_temporal = false;
        M_infty_ = CONSTRUCT_M_INFTY_(s, norms, T, w_infty, scaling_non_temporal);
        C_no_exp_all = R_' * M_ * R_;
        D = R_' * M_infty_ * R_;

        % Pressure component
        norms = [1, 0, 0, 0];
        M_ = CONSTRUCT_M_(s, norms, w_infty);
        C_no_exp_p = R_' * M_ * R_;

        % Kinetic component
        norms = [0, 1, 1, 0];
        M_ = CONSTRUCT_M_(s, norms, w_infty);
        C_no_exp_k = R_' * M_ * R_;

        % Entropy component
        norms = [0, 0, 0, 1];
        M_ = CONSTRUCT_M_(s, norms, w_infty);
        C_no_exp_S = R_' * M_ * R_;

    else
        R = CONSTRUCT_R(s);
        norms = [1, 1, 1, 1];
        M = CONSTRUCT_M(s, norms);
        M_hat = CONSTRUCT_M_HAT(s, norms);
        C_no_exp_all = R' * M * R;
        D = R' * M_hat * R;

        % Pressure component
        norms = [1, 0, 0, 0];
        M = CONSTRUCT_M(s, norms);
        C_no_exp_p = R' * M * R;

        % Kinetic component
        norms = [0, 1, 1, 0];
        M = CONSTRUCT_M(s, norms);
        C_no_exp_k = R' * M * R;

        % Entropy component
        norms = [0, 0, 0, 1];
        M = CONSTRUCT_M(s, norms);
        C_no_exp_S = R' * M * R;
    end

    %% Detect GPU availability and pre-allocate
    use_gpu = (gpuDeviceCount > 0);
    n_vec = size(A_, 1);

    if use_gpu
        gpu = gpuDevice;
        disp(gpu.Name + " selected.")
        gA_            = gpuArray(A_ * dt);
        gC_no_exp_all  = gpuArray(C_no_exp_all);
        gC_no_exp_p    = gpuArray(C_no_exp_p);
        gC_no_exp_k    = gpuArray(C_no_exp_k);
        gC_no_exp_S    = gpuArray(C_no_exp_S);
        v_buf          = gpuArray(zeros(n_vec, 1));
        vout_buf       = gpuArray(zeros(n_vec, 1));
        term_buf       = gpuArray(zeros(n_vec, 1));
    else
        disp("No CUDA device found. Running on CPU.")
        gA_            = A_ * dt;
        gC_no_exp_all  = C_no_exp_all;
        gC_no_exp_p    = C_no_exp_p;
        gC_no_exp_k    = C_no_exp_k;
        gC_no_exp_S    = C_no_exp_S;
        v_buf          = zeros(n_vec, 1);
        vout_buf       = zeros(n_vec, 1);
        term_buf       = zeros(n_vec, 1);
    end

    %% Integration (GPU or CPU)
    if get_amplification_only
        tic
        [t, E_all, E_p, E_k, E_S] = LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(v, n_t, dt, 0, orderTaylor, gA_, gC_no_exp_all, gC_no_exp_p, gC_no_exp_k, gC_no_exp_S, use_gpu, v_buf, vout_buf, term_buf);
        toc
    else
        tic
        [t, E_all, E_p, E_k, E_S, vf] = LINEAR_INTEGRATION_AND_ENERGY_GPU(v, n_t, N_samples, N_skip, dt, 0, orderTaylor, gA_, gC_no_exp_all, gC_no_exp_p, gC_no_exp_k, gC_no_exp_S, use_gpu, v_buf, vout_buf, term_buf);
        toc
    end

    %% GPU synchronization and memory cleanup
    if use_gpu
        wait(gpu);
        clear gA_ gC_no_exp_all gC_no_exp_p gC_no_exp_k gC_no_exp_S v_buf vout_buf term_buf;
        wait(gpu);
    else
        clear gA_ gC_no_exp_all gC_no_exp_p gC_no_exp_k gC_no_exp_S v_buf vout_buf term_buf;
    end

    %% Set initial energy values
    E_all(1) = real(v)' * C_no_exp_all * real(v);
    E_p(1) = real(v)' * C_no_exp_p * real(v);
    E_k(1) = real(v)' * C_no_exp_k * real(v);
    E_S(1) = real(v)' * C_no_exp_S * real(v);

    %% Plotting parameters
    fontsize_labels = 16;
    fontsize_ticks = 16;
    fontsize_legend = 16;
    interval_markers = 20;

    %% Energy gain plots
    if freestream_disturbances
        dE_dt_inflow_time = v' * D * v;
        T_ref = GET_T_REF(s);
        scaling_E_ref = real(dE_dt_inflow_time * T_ref);
        scaling_E = real(dE_dt_inflow_time .* t);
    else
        scaling_E_ref = real(E_all(1));
        scaling_E = real(E_all(1));
    end

    % Reference-scaled energy gain figure (freestream only)
    figure()
    hold on
    plot(t, real(E_all) ./ scaling_E_ref, 'r-', 'LineWidth', 1.5, 'DisplayName', '$E$');
    plot(t, real(E_p) ./ scaling_E_ref, 'g--', 'LineWidth', 1.5, 'DisplayName', '$E_p$');
    plot(t, real(E_k) ./ scaling_E_ref, 'b--', 'LineWidth', 1.5, 'DisplayName', '$E_k$');
    plot(t, real(E_S) ./ scaling_E_ref, 'c--', 'LineWidth', 1.5, 'DisplayName', '$E_s$');
    legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'best')
    xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
    if freestream_disturbances
        ylabel('$\displaystyle \frac{E(t)}{E_\infty^{ref}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
    else
        ylabel('$\displaystyle \frac{E(t)}{E(0)}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
    end
    box on
    xlim([0, T_f])
    set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')

    if ~get_amplification_only
        %% Energy time derivatives
        dE_all = real(E_all(3:end, 1) - E_all(1:end-2, 1)) / dt / 2;
        dE_p = real(E_p(3:end, 1) - E_p(1:end-2, 1)) / dt / 2;
        dE_k = real(E_k(3:end, 1) - E_k(1:end-2, 1)) / dt / 2;
        dE_S = real(E_S(3:end, 1) - E_S(1:end-2, 1)) / dt / 2;
        t_sample = t(N_skip+1:N_skip:end, 1);

        % Ensure consistent sampling for all arrays
        sample_dE_indices = N_skip:N_skip:length(dE_k);
        sample_dE_indices = sample_dE_indices(1:min(N_samples, length(sample_dE_indices)));
        t_sample = t_sample(1:length(sample_dE_indices));
        N_samples_actual = length(sample_dE_indices);

        if freestream_disturbances
            [dE_dt_inflow, ~] = ENERGY_INFLOW(v, s, w_infty);
            scaling_dE_dt = abs(dE_dt_inflow);
        else
            scaling_dE_dt = abs(dE_all(1, 1));
        end

        %% Compute energy budgets
        A_adv = zeros(N_samples_actual, 1);
        P_mom = zeros(N_samples_actual, 1);
        P_mass = zeros(N_samples_actual, 1);
        Dilat_P = zeros(N_samples_actual, 1);
        Transport_p = zeros(N_samples_actual, 1);
        Transport_tau = zeros(N_samples_actual, 1);
        Dissip = zeros(N_samples_actual, 1);
        A_adv_S = zeros(N_samples_actual, 1);
        P_mom_S = zeros(N_samples_actual, 1);
        P_mass_S = zeros(N_samples_actual, 1);
        P_T_S = zeros(N_samples_actual, 1);
        Transport_S = zeros(N_samples_actual, 1);
        Dissip_S = zeros(N_samples_actual, 1);
        Source_S = zeros(N_samples_actual, 1);
        shock_adv_acoustic = zeros(N_samples_actual, 1);
        shock_adv_kinetic = zeros(N_samples_actual, 1);
        shock_work_kinetic = zeros(N_samples_actual, 1);
        shock_adv_entropic = zeros(N_samples_actual, 1);

        output_flow = true;

        for i = 1:N_samples_actual
            budgets_kinetic = COMPUTE_BUDGETS_KINETIC(vf(:, i), s, output_flow);
            budgets_entropy = COMPUTE_BUDGETS_ENTROPY(vf(:, i), s, chemistry, output_flow);
            budgets_shock = COMPUTE_BUDGETS_SHOCK(vf(:, i), w_infty, s, chemistry);

            A_adv(i, 1) = budgets_kinetic.A_adv_sum;
            P_mom(i, 1) = budgets_kinetic.P_mom_sum;
            P_mass(i, 1) = budgets_kinetic.P_mass_sum;
            Dilat_P(i, 1) = budgets_kinetic.Dilat_P_sum;
            Dissip(i, 1) = budgets_kinetic.Dissipation_sum;
            Transport_p(i, 1) = budgets_kinetic.Transport_p_sum;
            Transport_tau(i, 1) = budgets_kinetic.Transport_tau_sum;
            A_adv_S(i, 1) = budgets_entropy.A_adv_sum;
            P_mom_S(i, 1) = budgets_entropy.P_mom_sum;
            P_mass_S(i, 1) = budgets_entropy.P_mass_sum;
            P_T_S(i, 1) = budgets_entropy.P_T_sum;
            Transport_S(i, 1) = budgets_entropy.Transport_sum;
            Dissip_S(i, 1) = budgets_entropy.Dissipation_sum;
            Source_S(i, 1) = budgets_entropy.Source_sum;
            shock_adv_acoustic(i, 1) = budgets_shock.adv_acoustic;
            shock_adv_kinetic(i, 1) = budgets_shock.adv_kinetic;
            shock_work_kinetic(i, 1) = budgets_shock.work_kinetic;
            shock_adv_entropic(i, 1) = budgets_shock.adv_entropic;
        end

        % %% Plot energy rate and shock budget terms
        % figure()
        % plot(t(2:end-1, 1), dE_all / scaling_dE_dt, 'r-', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}$');
        % hold on
        % plot(t(2:end-1, 1), dE_p / scaling_dE_dt, 'g--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}_p$');
        % plot(t(2:end-1, 1), dE_k / scaling_dE_dt, 'b--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}_k$');
        % plot(t(2:end-1, 1), dE_S / scaling_dE_dt, 'c--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}_s$');
        % plot(t_sample, shock_adv_acoustic / scaling_dE_dt, 'g-', 'LineWidth', 1.5, 'DisplayName', '$Shock-acoustic$');
        % plot(t_sample, shock_adv_kinetic / scaling_dE_dt, 'b-', 'LineWidth', 1.5, 'DisplayName', '$Shock-kinetic$');
        % plot(t_sample, shock_work_kinetic / scaling_dE_dt, 'b-*', 'LineWidth', 1.5, 'DisplayName', '$Shock-work-kinetic$');
        % plot(t_sample, shock_adv_entropic / scaling_dE_dt, 'c-', 'LineWidth', 1.5, 'DisplayName', '$Shock-entropic$');
        % yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        % legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'best')
        % xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % if freestream_disturbances
        %     ylabel('$\displaystyle \frac{\dot{E}(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % else
        %     ylabel('$\displaystyle \frac{\dot{E}(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % end
        % box on
        % set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        % hold off
        % 
        % %% Kinetic energy budget -- diagnostic figure
        % figure()
        % hold on
        % plot(t(2:end-1, 1), dE_k / scaling_dE_dt, 'b--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}^k$');
        % plot(t_sample, P_mom / scaling_dE_dt, 'm--', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{P}_u^k$');
        % plot(t_sample, P_mass / scaling_dE_dt, 'm-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{P}_\rho^k$');
        % plot(t_sample, Dilat_P / scaling_dE_dt, 'k-', 'LineWidth', 1.5, 'DisplayName', '$\Pi_d^k$');
        % plot(t_sample, (Transport_p + Transport_tau) / scaling_dE_dt, 'g-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{T}^k$');
        % plot(t_sample, A_adv / scaling_dE_dt, 'r--', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{A}^k$');
        % plot(t_sample, Dissip / scaling_dE_dt, 'c-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{D}_u^k$');
        % plot(t_sample, (Dissip + Dilat_P + P_mom + P_mass + Transport_p + Transport_tau + A_adv) / scaling_dE_dt, 'g--', 'LineWidth', 1.5, 'DisplayName', '$SUM^k$');
        % plot(t_sample, shock_adv_kinetic / scaling_dE_dt, 'b-', 'LineWidth', 1.5, 'DisplayName', '$Shock-kinetic$');
        % plot(t_sample, shock_work_kinetic / scaling_dE_dt, 'b-*', 'LineWidth', 1.5, 'DisplayName', '$Shock-work-kinetic$');
        % legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'best')
        % xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % if freestream_disturbances
        %     ylabel('$\displaystyle \frac{\dot{E}^k(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % else
        %     ylabel('$\displaystyle \frac{\dot{E}^k(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % end
        % box on
        % set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        % hold off

        %% Kinetic energy budget -- publication figure (grouped)
        figure()
        hold on
        plot(t(2:end-1, 1), dE_k / scaling_dE_dt, 'b--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}^k$');
        plot(t_sample, (P_mom + P_mass) / scaling_dE_dt, 'm--o', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{P}^k$');
        plot(t_sample, Dilat_P / scaling_dE_dt, 'r--s', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\Pi_d^k$');
        plot(t_sample, (A_adv + Transport_p + Transport_tau) / scaling_dE_dt, ...
            '-.^', 'Color', [1 0.5 0], 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{A}^k + \mathcal{T}^k$');
        plot(t_sample, Dissip / scaling_dE_dt, 'k:v', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{D}^k$');
        yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'best')
        xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        if freestream_disturbances
            ylabel('$\displaystyle \frac{\dot{E}^k(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        else
            ylabel('$\displaystyle \frac{\dot{E}^k(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        end
        box on
        set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        hold off

        %% Kinetic energy budget -- publication figure (detailed)
        figure()
        hold on
        plot(t(2:end-1, 1), dE_k / scaling_dE_dt, 'b--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}^k$');
        plot(t_sample, P_mom / scaling_dE_dt, 'm:d', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{P}_u^k$');
        plot(t_sample, Transport_p / scaling_dE_dt, ...
            '-o', 'Color', [1 0.5 0], 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{T}^k_p$');
        plot(t_sample, Dissip / scaling_dE_dt, 'k:v', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{D}^k$');
        plot(t_sample, A_adv / scaling_dE_dt, '-.d', 'Color', [0.75 0.75 0], 'LineWidth', 1.5, ...
            'MarkerIndices', interval_markers/2:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{A}^k$');
        yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'best')
        xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        if freestream_disturbances
            ylabel('$\displaystyle \frac{\dot{E}^k(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        else
            ylabel('$\displaystyle \frac{\dot{E}^k(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        end
        box on
        set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        hold off

        % %% Entropic energy budget -- diagnostic figure
        % figure()
        % hold on
        % plot(t(2:end-1, 1), dE_S / scaling_dE_dt, 'c--', 'LineWidth', 1.5, 'DisplayName', 'Total entropic');
        % plot(t_sample, P_mom_S / scaling_dE_dt, 'm--', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{P}_u^s$');
        % plot(t_sample, P_mass_S / scaling_dE_dt, 'm:d', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{P}_\rho^s$');
        % plot(t_sample, P_T_S / scaling_dE_dt, 'k-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{P}_T^s$');
        % plot(t_sample, (Transport_S + A_adv_S) / scaling_dE_dt, 'g-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{A}^s + \mathcal{T}^s$');
        % plot(t_sample, Dissip_S / scaling_dE_dt, 'c-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{D}^s$');
        % plot(t_sample, Source_S / scaling_dE_dt, 'r-', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{S}^s$');
        % plot(t_sample, (Dissip_S + P_mom_S + P_mass_S + P_T_S + Transport_S + A_adv_S + Source_S) / scaling_dE_dt, 'g--', 'LineWidth', 1.5, 'DisplayName', '$SUM^s$');
        % plot(t_sample, shock_adv_entropic / scaling_dE_dt, 'b-', 'LineWidth', 1.5, 'DisplayName', '$Shock-entropic$');
        % legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'southwest')
        % xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % if freestream_disturbances
        %     ylabel('$\displaystyle \frac{\dot{E}^s(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % else
        %     ylabel('$\displaystyle \frac{\dot{E}^s(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        % end
        % box on
        % set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        % hold off

        %% Entropic energy budget -- publication figure (grouped)
        figure()
        hold on
        plot(t(2:end-1, 1), dE_S / scaling_dE_dt, 'c--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}^s$');
        plot(t_sample, (P_mom_S + P_mass_S + P_T_S) / scaling_dE_dt, 'm--o', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{P}^s$');
        plot(t_sample, (Transport_S + A_adv_S) / scaling_dE_dt, '-.^', 'Color', [1 0.5 0], 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{A}^s + \mathcal{T}^s$');
        plot(t_sample, (Dissip_S) / scaling_dE_dt, 'k:v', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{D}^s$');
        plot(t_sample, Source_S / scaling_dE_dt, 'r-p', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{S}^s$');
        yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'southwest')
        xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        if freestream_disturbances
            ylabel('$\displaystyle \frac{\dot{E}^s(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        else
            ylabel('$\displaystyle \frac{\dot{E}^s(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        end
        box on
        set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        hold off

        %% Entropic energy budget -- publication figure (detailed)
        figure()
        hold on
        plot(t(2:end-1, 1), dE_S / scaling_dE_dt, 'c--', 'LineWidth', 1.5, 'DisplayName', '$\dot{E}^s$');
        plot(t_sample, P_mom_S / scaling_dE_dt, 'm:d', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{P}_u^s$');
        plot(t_sample, (Dissip_S) / scaling_dE_dt, ...
            'k:v', 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{D}^s$');
        plot(t_sample, A_adv_S / scaling_dE_dt, '-.d', 'Color', [1 0.5 0], 'LineWidth', 1.5, ...
            'MarkerIndices', 1:interval_markers:length(t_sample), 'MarkerSize', 5, 'DisplayName', '$\mathcal{A}^s$');
        yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
        legend('Interpreter', 'latex', 'FontSize', fontsize_legend, 'Location', 'southwest')
        xlabel('$t U_\infty / R$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        if freestream_disturbances
            ylabel('$\displaystyle \frac{\dot{E}^s(t)}{\overline{\dot{E}_\infty}}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        else
            ylabel('$\displaystyle \frac{\dot{E}^s(t)}{|\dot{E}(0)|}$', 'Interpreter', 'latex', 'FontSize', fontsize_labels)
        end
        box on
        set(gca, 'FontSize', fontsize_ticks, 'TickLabelInterpreter', 'latex')
        hold off

        %% Store energetic budget gains
        max_gain.budgets.max_prod_u = max(P_mom / scaling_dE_dt, [], "all");
        max_gain.budgets.max_transport = max(abs(Transport_p) / scaling_dE_dt, [], "all");
        max_gain.budgets.start_adv = A_adv(1, 1) / scaling_dE_dt;
        max_gain.budgets.start_transport = Transport_p(1, 1) / scaling_dE_dt;

        max_gain.budgets.dE_p_ref = mean(real(shock_adv_acoustic / scaling_dE_dt));
        max_gain.budgets.dE_k_ref = mean(real((shock_adv_kinetic + shock_work_kinetic) / scaling_dE_dt));
        max_gain.budgets.dE_s_ref = mean(real(shock_adv_entropic / scaling_dE_dt));
    end

    %% Compute maximum amplification factors
    start_index = 10;
    max_gain.temporal.all = max(real(E_all(start_index:end) ./ scaling_E(start_index:end)), [], "all");
    max_gain.temporal.acoustic = max(real(E_p(start_index:end) ./ scaling_E(start_index:end)), [], "all");
    max_gain.temporal.entropic = max(real(E_S(start_index:end) ./ scaling_E(start_index:end)), [], "all");
    max_gain.temporal.kinetic = max(real(E_k(start_index:end) ./ scaling_E(start_index:end)), [], "all");

    max_gain.non_temporal.all = max(real(E_all ./ scaling_E_ref), [], "all");
    max_gain.non_temporal.acoustic = max(real(E_p ./ scaling_E_ref), [], "all");
    max_gain.non_temporal.entropic = max(real(E_S ./ scaling_E_ref), [], "all");
    max_gain.non_temporal.kinetic = max(real(E_k ./ scaling_E_ref), [], "all");
end


function [t, E_all, E_p, E_k, E_S, Vf] = LINEAR_INTEGRATION_AND_ENERGY_GPU(x, n_t, N_samples, N_skip, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)
% LINEAR_INTEGRATION_AND_ENERGY_GPU  Forward integration with energy tracking and sampling.
%   Works on both GPU and CPU using pre-allocated buffers.
%
%   [t, E_all, E_p, E_k, E_S, Vf] = LINEAR_INTEGRATION_AND_ENERGY_GPU(...)
%
%   Inputs:
%       x             - Initial state vector (CPU)
%       n_t           - Number of time steps
%       N_samples     - Number of state snapshots to store
%       N_skip        - Steps between successive snapshots
%       dt            - Time step size
%       t_0           - Initial time
%       orderTaylor   - Taylor series order for matrix exponential
%       Adt           - System matrix A*dt (gpuArray or CPU sparse)
%       C_no_exp_all  - Total energy weight operator (gpuArray or CPU sparse)
%       C_no_exp_p    - Pressure energy weight operator
%       C_no_exp_k    - Kinetic energy weight operator
%       C_no_exp_S    - Entropy energy weight operator
%       use_gpu       - Logical flag: true for GPU path, false for CPU path
%       v_buf         - Pre-allocated buffer for working vector
%       vout_buf      - Pre-allocated buffer for output accumulation
%       term_buf      - Pre-allocated buffer for Taylor term
%
%   Outputs:
%       t     - Time vector (CPU)
%       E_all - Total energy history (CPU)
%       E_p   - Pressure energy history (CPU)
%       E_k   - Kinetic energy history (CPU)
%       E_S   - Entropy energy history (CPU)
%       Vf    - Sampled state vectors (CPU, columns)

    %% Allocate storage
    E_all = zeros(n_t + 1, 1);
    E_p = zeros(n_t + 1, 1);
    E_k = zeros(n_t + 1, 1);
    E_S = zeros(n_t + 1, 1);
    t = zeros(n_t + 1, 1);
    n_vec = numel(x);
    Vf_store = zeros(n_vec, N_samples);

    t(1) = t_0;

    %% Copy initial vector into pre-allocated buffer
    v_buf(:) = x;
    vout_buf(:) = v_buf;

    %% Forward propagation via Taylor-series exponential (no allocation)
    count = 1;
    for i = 1:n_t
        factorial_val = 1;
        for j = 1:orderTaylor
            term_buf(:) = Adt * v_buf;
            v_buf(:) = term_buf;
            factorial_val = j * factorial_val;
            vout_buf(:) = vout_buf + v_buf / factorial_val;
        end
        v_buf(:) = vout_buf;
        if mod(i, N_skip) == 0 && count <= N_samples
            if use_gpu
                Vf_store(:, count) = gather(v_buf);
            else
                Vf_store(:, count) = v_buf;
            end
            count = count + 1;
        end

        % Store iteration data
        t(i+1) = t(i) + dt;
        if use_gpu
            v_real = gather(real(v_buf));
        else
            v_real = real(v_buf);
        end
        E_all(i+1) = v_real' * gather(C_no_exp_all * real(v_buf));
        E_p(i+1) = v_real' * gather(C_no_exp_p * real(v_buf));
        E_k(i+1) = v_real' * gather(C_no_exp_k * real(v_buf));
        E_S(i+1) = v_real' * gather(C_no_exp_S * real(v_buf));
    end

    Vf = Vf_store;
end


function [t, E_all, E_p, E_k, E_S] = LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(x, n_t, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)
% LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES  Forward integration (energy only, no snapshots).
%   Works on both GPU and CPU using pre-allocated buffers.
%
%   [t, E_all, E_p, E_k, E_S] = LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(...)
%
%   Same as LINEAR_INTEGRATION_AND_ENERGY_GPU but without storing intermediate
%   state snapshots, reducing memory usage.
%
%   Inputs:
%       x             - Initial state vector (CPU)
%       n_t           - Number of time steps
%       dt            - Time step size
%       t_0           - Initial time
%       orderTaylor   - Taylor series order for matrix exponential
%       Adt           - System matrix A*dt (gpuArray or CPU sparse)
%       C_no_exp_all  - Total energy weight operator
%       C_no_exp_p    - Pressure energy weight operator
%       C_no_exp_k    - Kinetic energy weight operator
%       C_no_exp_S    - Entropy energy weight operator
%       use_gpu       - Logical flag: true for GPU path, false for CPU path
%       v_buf         - Pre-allocated buffer for working vector
%       vout_buf      - Pre-allocated buffer for output accumulation
%       term_buf      - Pre-allocated buffer for Taylor term
%
%   Outputs:
%       t     - Time vector (CPU)
%       E_all - Total energy history (CPU)
%       E_p   - Pressure energy history (CPU)
%       E_k   - Kinetic energy history (CPU)
%       E_S   - Entropy energy history (CPU)

    %% Allocate storage
    E_all = zeros(n_t + 1, 1);
    E_p = zeros(n_t + 1, 1);
    E_k = zeros(n_t + 1, 1);
    E_S = zeros(n_t + 1, 1);
    t = zeros(n_t + 1, 1);

    t(1) = t_0;

    %% Copy initial vector into pre-allocated buffer
    v_buf(:) = x;
    vout_buf(:) = v_buf;

    %% Forward propagation via Taylor-series exponential (no allocation)
    for i = 1:n_t
        factorial_val = 1;
        for j = 1:orderTaylor
            term_buf(:) = Adt * v_buf;
            v_buf(:) = term_buf;
            factorial_val = j * factorial_val;
            vout_buf(:) = vout_buf + v_buf / factorial_val;
        end
        v_buf(:) = vout_buf;

        % Store iteration data
        t(i+1) = t(i) + dt;
        if use_gpu
            v_real = gather(real(v_buf));
        else
            v_real = real(v_buf);
        end
        E_all(i+1) = v_real' * gather(C_no_exp_all * real(v_buf));
        E_p(i+1) = v_real' * gather(C_no_exp_p * real(v_buf));
        E_k(i+1) = v_real' * gather(C_no_exp_k * real(v_buf));
        E_S(i+1) = v_real' * gather(C_no_exp_S * real(v_buf));
    end
end
