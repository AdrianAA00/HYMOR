using LinearAlgebra
using SparseArrays
using CUDA
using CUDA.CUSPARSE
using Printf
using Plots

"""
    LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances, v, lambda, T_opt, A_, s, chemistry, T_f, get_amplification_only, w_infty)

Integrate the linearized system, plot energy gains, and compute budgets.

Propagates the optimal initial perturbation v forward in time under the
linearized dynamics defined by A_ using a GPU-accelerated Taylor-series
matrix exponential. Computes the decomposed energy evolution (total,
pressure, kinetic, entropy) and generates publication-quality figures
of energy gains and budgets.

# Arguments
- `freestream_disturbances`: If true, use the extended system (downstream + freestream).
- `v`: Optimal initial perturbation eigenvector.
- `lambda`: Optimal gain (eigenvalue) for overlay on plots.
- `T_opt`: Optimal time horizon corresponding to lambda.
- `A_`: System matrix (or extended system matrix).
- `s`: Solution Dict{String,Any} with flow fields and mesh data.
- `chemistry`: Chemistry model Dict{String,Any}.
- `T_f`: Final integration time for plotting.
- `get_amplification_only`: If true, skip budget computations and sampling.
- `w_infty`: Vector of freestream disturbance frequencies.

# Returns
- `max_gain`: Dict{String,Any} with fields:
    - `"temporal"` => Dict with `"all"`, `"acoustic"`, `"kinetic"`, `"entropic"`
    - `"non_temporal"` => Dict with `"all"`, `"acoustic"`, `"kinetic"`, `"entropic"`
    - `"budgets"` (only when !get_amplification_only)

# Notes
- Uses 5th-order Taylor expansion for the matrix exponential.
- All heavy computations are performed on the GPU.
- Generates multiple publication-quality figures for energy gains,
  kinetic energy budgets, and entropic energy budgets.

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function LINEAR_INTEGRATION_AND_GAINS(freestream_disturbances::Bool, v, lambda, T_opt, A_, s::Dict{String,Any},
                                       chemistry::Dict{String,Any}, T_f, get_amplification_only::Bool, w_infty)
    SET_PLOT_DEFAULTS()

    ## Time-stepping parameters
    orderTaylor = 5  # Order of Taylor expansion for matrix exponential approximation

    s = CFL_TIMESTEP(s)
    t_temp = s["time_integration"]["dt"]
    n_t = round(Int, T_f / t_temp)
    dt = T_f / n_t
    N_skip = div(s["mesh"]["Neta"], 2)
    N_samples = max(div(n_t, N_skip), 500)

    println()
    @printf("Number timesteps = %f \n", n_t)

    ## Construct energy-norm operators
    if freestream_disturbances
        R_ = CONSTRUCT_R_(s, w_infty)
        norms = [1, 1, 1, 1]
        M_ = CONSTRUCT_M_(s, norms, w_infty)

        T_val = 1
        scaling_non_temporal = false
        M_infty_ = CONSTRUCT_M_INFTY_(s, norms, T_val, w_infty, scaling_non_temporal)
        C_no_exp_all = R_' * M_ * R_
        D_mat = R_' * M_infty_ * R_

        # Pressure component
        norms = [1, 0, 0, 0]
        M_ = CONSTRUCT_M_(s, norms, w_infty)
        C_no_exp_p = R_' * M_ * R_

        # Kinetic component
        norms = [0, 1, 1, 0]
        M_ = CONSTRUCT_M_(s, norms, w_infty)
        C_no_exp_k = R_' * M_ * R_

        # Entropy component
        norms = [0, 0, 0, 1]
        M_ = CONSTRUCT_M_(s, norms, w_infty)
        C_no_exp_S = R_' * M_ * R_

    else
        R = CONSTRUCT_R(s)
        norms = [1, 1, 1, 1]
        M = CONSTRUCT_M(s, norms)
        M_hat = CONSTRUCT_M_HAT(s, norms)
        C_no_exp_all = R' * M * R
        D_mat = R' * M_hat * R

        # Pressure component
        norms = [1, 0, 0, 0]
        M = CONSTRUCT_M(s, norms)
        C_no_exp_p = R' * M * R

        # Kinetic component
        norms = [0, 1, 1, 0]
        M = CONSTRUCT_M(s, norms)
        C_no_exp_k = R' * M * R

        # Entropy component
        norms = [0, 0, 0, 1]
        M = CONSTRUCT_M(s, norms)
        C_no_exp_S = R' * M * R
    end

    ## Detect GPU availability and pre-allocate
    use_gpu = CUDA.functional()
    n_vec = size(A_, 1)

    if use_gpu
        gpu = CUDA.device()
        println("$(CUDA.name(gpu)) selected.")
        gA_           = CuSparseMatrixCSC(A_ * dt)
        gC_no_exp_all = CuSparseMatrixCSC(C_no_exp_all)
        gC_no_exp_p   = CuSparseMatrixCSC(C_no_exp_p)
        gC_no_exp_k   = CuSparseMatrixCSC(C_no_exp_k)
        gC_no_exp_S   = CuSparseMatrixCSC(C_no_exp_S)
        v_buf         = CUDA.zeros(ComplexF64, n_vec) # Can be complex in freestream receptivity case
        vout_buf      = CUDA.zeros(ComplexF64, n_vec)
        term_buf      = CUDA.zeros(ComplexF64, n_vec)
    else
        println("No CUDA device found. Running on CPU.")
        gA_           = A_ * dt
        gC_no_exp_all = C_no_exp_all
        gC_no_exp_p   = C_no_exp_p
        gC_no_exp_k   = C_no_exp_k
        gC_no_exp_S   = C_no_exp_S
        v_buf         = zeros(ComplexF64, n_vec)
        vout_buf      = zeros(ComplexF64, n_vec)
        term_buf      = zeros(ComplexF64, n_vec)
    end

    ## Integration (GPU or CPU)
    if get_amplification_only
        time_start = time()
        (t, E_all, E_p, E_k, E_S) = LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(v, n_t, dt, 0.0, orderTaylor, gA_, gC_no_exp_all, gC_no_exp_p, gC_no_exp_k, gC_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)
        elapsed = time() - time_start
        @printf("Elapsed time: %.6f seconds\n", elapsed)
    else
        time_start = time()
        (t, E_all, E_p, E_k, E_S, vf) = LINEAR_INTEGRATION_AND_ENERGY_GPU(v, n_t, N_samples, N_skip, dt, 0.0, orderTaylor, gA_, gC_no_exp_all, gC_no_exp_p, gC_no_exp_k, gC_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)
        elapsed = time() - time_start
        @printf("Elapsed time: %.6f seconds\n", elapsed)
    end

    ## GPU synchronization and memory cleanup
    if use_gpu
        CUDA.synchronize()
        CUDA.unsafe_free!(gA_)
        CUDA.unsafe_free!(gC_no_exp_all)
        CUDA.unsafe_free!(gC_no_exp_p)
        CUDA.unsafe_free!(gC_no_exp_k)
        CUDA.unsafe_free!(gC_no_exp_S)
        CUDA.unsafe_free!(v_buf)
        CUDA.unsafe_free!(vout_buf)
        CUDA.unsafe_free!(term_buf)
        CUDA.synchronize()
        CUDA.reclaim()
    end

    ## Set initial energy values
    v_real = real.(v)
    E_all[1] = real(v_real' * C_no_exp_all * v_real)
    E_p[1] = real(v_real' * C_no_exp_p * v_real)
    E_k[1] = real(v_real' * C_no_exp_k * v_real)
    E_S[1] = real(v_real' * C_no_exp_S * v_real)

    ## Plotting parameters (font sizes handled by SET_PLOT_DEFAULTS)
    interval_markers = 20

    ## Energy gain plots
    if freestream_disturbances
        dE_dt_inflow_time = real(v' * D_mat * v)
        T_ref = GET_T_REF(s)
        scaling_E_ref = real(dE_dt_inflow_time * T_ref)
        scaling_E = real.(dE_dt_inflow_time .* t)
    else
        scaling_E_ref = real(E_all[1])
        scaling_E = real(E_all[1])
    end

    # Reference-scaled energy gain figure
    p1 = plot(t, real.(E_all) ./ scaling_E_ref, color=:red, linestyle=:solid, linewidth=1.5, label="\$E\$")
    plot!(p1, t, real.(E_p) ./ scaling_E_ref, color=:green, linestyle=:dash, linewidth=1.5, label="\$E_p\$")
    plot!(p1, t, real.(E_k) ./ scaling_E_ref, color=:blue, linestyle=:dash, linewidth=1.5, label="\$E_k\$")
    plot!(p1, t, real.(E_S) ./ scaling_E_ref, color=:cyan, linestyle=:dash, linewidth=1.5, label="\$E_s\$")
    if freestream_disturbances
        ylabel!(p1, "\$E(t)/E_{\\infty}^{ref}\$")
    else
        ylabel!(p1, "\$E(t)/E(0)\$")
    end
    xlabel!(p1, "\$t U_{\\infty} / R\$")
    xlims!(p1, (0, T_f))
    mkpath(_FIG_OUTPUT_DIR[])
    Plots.savefig(p1, joinpath(_FIG_OUTPUT_DIR[], "gains_energy_out.png"))
    display(p1)

    ## Initialize max_gain dict
    max_gain = Dict{String,Any}()
    max_gain["temporal"] = Dict{String,Any}()
    max_gain["non_temporal"] = Dict{String,Any}()

    if !get_amplification_only
        ## Energy time derivatives
        dE_all = real.(E_all[3:end] .- E_all[1:end-2]) ./ dt ./ 2
        dE_p = real.(E_p[3:end] .- E_p[1:end-2]) ./ dt ./ 2
        dE_k = real.(E_k[3:end] .- E_k[1:end-2]) ./ dt ./ 2
        dE_S = real.(E_S[3:end] .- E_S[1:end-2]) ./ dt ./ 2
        t_sample = t[N_skip+1:N_skip:end]

        # Ensure consistent sampling for all arrays
        sample_dE_indices = collect(N_skip:N_skip:length(dE_k))
        sample_dE_indices = sample_dE_indices[1:min(N_samples, length(sample_dE_indices))]
        t_sample = t_sample[1:length(sample_dE_indices)]
        N_samples_actual = length(sample_dE_indices)

        if freestream_disturbances
            (dE_dt_inflow, _) = ENERGY_INFLOW(v, s, w_infty)
            scaling_dE_dt = abs(dE_dt_inflow)
        else
            scaling_dE_dt = abs(dE_all[1])
        end

        ## Compute energy budgets
        A_adv = zeros(N_samples_actual, 1)
        P_mom = zeros(N_samples_actual, 1)
        P_mass = zeros(N_samples_actual, 1)
        Dilat_P = zeros(N_samples_actual, 1)
        Transport_p = zeros(N_samples_actual, 1)
        Transport_tau = zeros(N_samples_actual, 1)
        Dissip = zeros(N_samples_actual, 1)
        A_adv_S = zeros(N_samples_actual, 1)
        P_mom_S = zeros(N_samples_actual, 1)
        P_mass_S = zeros(N_samples_actual, 1)
        P_T_S = zeros(N_samples_actual, 1)
        Transport_S = zeros(N_samples_actual, 1)
        Dissip_S = zeros(N_samples_actual, 1)
        Source_S = zeros(N_samples_actual, 1)
        shock_adv_acoustic = zeros(N_samples_actual, 1)
        shock_adv_kinetic = zeros(N_samples_actual, 1)
        shock_work_kinetic = zeros(N_samples_actual, 1)
        shock_adv_entropic = zeros(N_samples_actual, 1)

        output_flow = true

        for i in 1:N_samples_actual
            budgets_kinetic = COMPUTE_BUDGETS_KINETIC(vf[:, i], s, output_flow)
            budgets_entropy = COMPUTE_BUDGETS_ENTROPY(vf[:, i], s, chemistry, output_flow)
            budgets_shock = COMPUTE_BUDGETS_SHOCK(vf[:, i], w_infty, s, chemistry)

            A_adv[i, 1] = budgets_kinetic["A_adv_sum"]
            P_mom[i, 1] = budgets_kinetic["P_mom_sum"]
            P_mass[i, 1] = budgets_kinetic["P_mass_sum"]
            Dilat_P[i, 1] = budgets_kinetic["Dilat_P_sum"]
            Dissip[i, 1] = budgets_kinetic["Dissipation_sum"]
            Transport_p[i, 1] = budgets_kinetic["Transport_p_sum"]
            Transport_tau[i, 1] = budgets_kinetic["Transport_tau_sum"]
            A_adv_S[i, 1] = budgets_entropy["A_adv_sum"]
            P_mom_S[i, 1] = budgets_entropy["P_mom_sum"]
            P_mass_S[i, 1] = budgets_entropy["P_mass_sum"]
            P_T_S[i, 1] = budgets_entropy["P_T_sum"]
            Transport_S[i, 1] = budgets_entropy["Transport_sum"]
            Dissip_S[i, 1] = budgets_entropy["Dissipation_sum"]
            Source_S[i, 1] = budgets_entropy["Source_sum"]
            shock_adv_acoustic[i, 1] = budgets_shock["adv_acoustic"]
            shock_adv_kinetic[i, 1] = budgets_shock["adv_kinetic"]
            shock_work_kinetic[i, 1] = budgets_shock["work_kinetic"]
            shock_adv_entropic[i, 1] = budgets_shock["adv_entropic"]
        end

        ## Kinetic energy budget -- publication figure (grouped)
        p2 = plot(t[2:end-1], dE_k ./ scaling_dE_dt, color=:blue, linestyle=:dash, linewidth=1.5, label="\$\\dot{E}^k\$")
        plot!(p2, t_sample, (P_mom .+ P_mass) ./ scaling_dE_dt, color=:magenta, linestyle=:dash, linewidth=1.5, markershape=:circle, label="\$\\mathcal{P}^k\$")
        plot!(p2, t_sample, Dilat_P ./ scaling_dE_dt, color=:red, linestyle=:dash, linewidth=1.5, markershape=:square, label="\$\\Pi_d^k\$")
        plot!(p2, t_sample, (A_adv .+ Transport_p .+ Transport_tau) ./ scaling_dE_dt, color=RGB(1, 0.5, 0), linestyle=:dashdot, linewidth=1.5, markershape=:utriangle, label="\$\\mathcal{A}^k + \\mathcal{T}^k\$")
        plot!(p2, t_sample, Dissip ./ scaling_dE_dt, color=:black, linestyle=:dot, linewidth=1.5, markershape=:dtriangle, label="\$\\mathcal{D}^k\$")
        hline!(p2, [0], color=:black, linewidth=0.5, label=nothing)
        if freestream_disturbances
            ylabel!(p2, "\$\\dot{E}^k(t)/\\overline{\\dot{E}_{\\infty}}\$")
        else
            ylabel!(p2, "\$\\dot{E}^k(t)/|\\dot{E}(0)|\$")
        end
        xlabel!(p2, "\$t U_{\\infty} / R\$")
        Plots.savefig(p2, joinpath(_FIG_OUTPUT_DIR[], "budget_kinetic_grouped_out.png"))
        display(p2)

        ## Kinetic energy budget -- publication figure (detailed)
        p3 = plot(t[2:end-1], dE_k ./ scaling_dE_dt, color=:blue, linestyle=:dash, linewidth=1.5, label="\$\\dot{E}^k\$")
        plot!(p3, t_sample, P_mom ./ scaling_dE_dt, color=:magenta, linestyle=:dot, linewidth=1.5, markershape=:diamond, label="\$\\mathcal{P}_u^k\$")
        plot!(p3, t_sample, Transport_p ./ scaling_dE_dt, color=RGB(1, 0.5, 0), linestyle=:solid, linewidth=1.5, markershape=:circle, label="\$\\mathcal{T}^k_p\$")
        plot!(p3, t_sample, Dissip ./ scaling_dE_dt, color=:black, linestyle=:dot, linewidth=1.5, markershape=:dtriangle, label="\$\\mathcal{D}^k\$")
        plot!(p3, t_sample, A_adv ./ scaling_dE_dt, color=RGB(0.75, 0.75, 0), linestyle=:dashdot, linewidth=1.5, markershape=:diamond, label="\$\\mathcal{A}^k\$")
        hline!(p3, [0], color=:black, linewidth=0.5, label=nothing)
        if freestream_disturbances
            ylabel!(p3, "\$\\dot{E}^k(t)/\\overline{\\dot{E}_{\\infty}}\$")
        else
            ylabel!(p3, "\$\\dot{E}^k(t)/|\\dot{E}(0)|\$")
        end
        xlabel!(p3, "\$t U_{\\infty} / R\$")
        Plots.savefig(p3, joinpath(_FIG_OUTPUT_DIR[], "budget_kinetic_detailed_out.png"))
        display(p3)

        ## Entropic energy budget -- publication figure (grouped)
        p4 = plot(t[2:end-1], dE_S ./ scaling_dE_dt, color=:cyan, linestyle=:dash, linewidth=1.5, label="\$\\dot{E}^s\$")
        plot!(p4, t_sample, (P_mom_S .+ P_mass_S .+ P_T_S) ./ scaling_dE_dt, color=:magenta, linestyle=:dash, linewidth=1.5, markershape=:circle, label="\$\\mathcal{P}^s\$")
        plot!(p4, t_sample, (Transport_S .+ A_adv_S) ./ scaling_dE_dt, color=RGB(1, 0.5, 0), linestyle=:dashdot, linewidth=1.5, markershape=:utriangle, label="\$\\mathcal{A}^s + \\mathcal{T}^s\$")
        plot!(p4, t_sample, Dissip_S ./ scaling_dE_dt, color=:black, linestyle=:dot, linewidth=1.5, markershape=:dtriangle, label="\$\\mathcal{D}^s\$")
        plot!(p4, t_sample, Source_S ./ scaling_dE_dt, color=:red, linestyle=:solid, linewidth=1.5, markershape=:pentagon, label="\$\\mathcal{S}^s\$")
        hline!(p4, [0], color=:black, linewidth=0.5, label=nothing)
        if freestream_disturbances
            ylabel!(p4, "\$\\dot{E}^s(t)/\\overline{\\dot{E}_{\\infty}}\$")
        else
            ylabel!(p4, "\$\\dot{E}^s(t)/|\\dot{E}(0)|\$")
        end
        xlabel!(p4, "\$t U_{\\infty} / R\$")
        Plots.savefig(p4, joinpath(_FIG_OUTPUT_DIR[], "budget_entropic_grouped_out.png"))
        display(p4)

        ## Entropic energy budget -- publication figure (detailed)
        p5 = plot(t[2:end-1], dE_S ./ scaling_dE_dt, color=:cyan, linestyle=:dash, linewidth=1.5, label="\$\\dot{E}^s\$")
        plot!(p5, t_sample, P_mom_S ./ scaling_dE_dt, color=:magenta, linestyle=:dot, linewidth=1.5, markershape=:diamond, label="\$\\mathcal{P}_u^s\$")
        plot!(p5, t_sample, Dissip_S ./ scaling_dE_dt, color=:black, linestyle=:dot, linewidth=1.5, markershape=:dtriangle, label="\$\\mathcal{D}^s\$")
        plot!(p5, t_sample, A_adv_S ./ scaling_dE_dt, color=RGB(1, 0.5, 0), linestyle=:dashdot, linewidth=1.5, markershape=:diamond, label="\$\\mathcal{A}^s\$")
        hline!(p5, [0], color=:black, linewidth=0.5, label=nothing)
        if freestream_disturbances
            ylabel!(p5, "\$\\dot{E}^s(t)/\\overline{\\dot{E}_{\\infty}}\$")
        else
            ylabel!(p5, "\$\\dot{E}^s(t)/|\\dot{E}(0)|\$")
        end
        xlabel!(p5, "\$t U_{\\infty} / R\$")
        Plots.savefig(p5, joinpath(_FIG_OUTPUT_DIR[], "budget_entropic_detailed_out.png"))
        display(p5)

        ## Store energetic budget gains
        max_gain["budgets"] = Dict{String,Any}()
        max_gain["budgets"]["max_prod_u"] = maximum(P_mom ./ scaling_dE_dt)
        max_gain["budgets"]["max_transport"] = maximum(abs.(Transport_p) ./ scaling_dE_dt)
        max_gain["budgets"]["start_adv"] = A_adv[1, 1] / scaling_dE_dt
        max_gain["budgets"]["start_transport"] = Transport_p[1, 1] / scaling_dE_dt

        max_gain["budgets"]["dE_p_ref"] = mean(real.(shock_adv_acoustic ./ scaling_dE_dt))
        max_gain["budgets"]["dE_k_ref"] = mean(real.((shock_adv_kinetic .+ shock_work_kinetic) ./ scaling_dE_dt))
        max_gain["budgets"]["dE_s_ref"] = mean(real.(shock_adv_entropic ./ scaling_dE_dt))
    end

    ## Compute maximum amplification factors
    start_index = 10
    max_gain["temporal"]["all"] = maximum(real.(E_all[start_index:end] ./ scaling_E[start_index:end]))
    max_gain["temporal"]["acoustic"] = maximum(real.(E_p[start_index:end] ./ scaling_E[start_index:end]))
    max_gain["temporal"]["entropic"] = maximum(real.(E_S[start_index:end] ./ scaling_E[start_index:end]))
    max_gain["temporal"]["kinetic"] = maximum(real.(E_k[start_index:end] ./ scaling_E[start_index:end]))

    max_gain["non_temporal"]["all"] = maximum(real.(E_all ./ scaling_E_ref))
    max_gain["non_temporal"]["acoustic"] = maximum(real.(E_p ./ scaling_E_ref))
    max_gain["non_temporal"]["entropic"] = maximum(real.(E_S ./ scaling_E_ref))
    max_gain["non_temporal"]["kinetic"] = maximum(real.(E_k ./ scaling_E_ref))

    return max_gain
end


"""
    LINEAR_INTEGRATION_AND_ENERGY_GPU(x, n_t, N_samples, N_skip, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)

Forward integration with energy tracking and sampling.
Works on both GPU and CPU using pre-allocated buffers.
"""
function LINEAR_INTEGRATION_AND_ENERGY_GPU(x, n_t, N_samples, N_skip, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu::Bool, v_buf, vout_buf, term_buf)
    ## Allocate storage (CPU)
    E_all = zeros(n_t + 1)
    E_p = zeros(n_t + 1)
    E_k = zeros(n_t + 1)
    E_S = zeros(n_t + 1)
    t = zeros(n_t + 1)
    n_elements = length(x)
    Vf = zeros(n_elements, N_samples)

    t[1] = t_0

    ## Copy initial vector into pre-allocated buffer
    copyto!(v_buf, x)
    copyto!(vout_buf, v_buf)

    ## Forward propagation via Taylor-series exponential (no allocation)
    count = 1
    for i in 1:n_t
        factorial_val = 1
        for j in 1:orderTaylor
            mul!(term_buf, Adt, v_buf, 1.0, 0.0)
            copyto!(v_buf, term_buf)
            factorial_val = j * factorial_val
            vout_buf .+= v_buf .* (1.0 / factorial_val)
        end
        copyto!(v_buf, vout_buf)
        if mod(i, N_skip) == 0 && count <= N_samples
            if use_gpu
                Vf[:, count] .= Array(v_buf)
            else
                Vf[:, count] .= v_buf
            end
            count += 1
        end

        # Store iteration data
        t[i+1] = t[i] + dt
        if use_gpu
            v_real = Array(real.(v_buf))
            v_real_dev = real.(v_buf)
        else
            v_real = real.(v_buf)
            v_real_dev = v_real
        end
        E_all[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_all * v_real_dev) : C_no_exp_all * v_real_dev))
        E_p[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_p * v_real_dev) : C_no_exp_p * v_real_dev))
        E_k[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_k * v_real_dev) : C_no_exp_k * v_real_dev))
        E_S[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_S * v_real_dev) : C_no_exp_S * v_real_dev))
    end

    return (t, E_all, E_p, E_k, E_S, Vf)
end


"""
    LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(x, n_t, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu, v_buf, vout_buf, term_buf)

Forward integration (energy only, no snapshots).
Works on both GPU and CPU using pre-allocated buffers.
"""
function LINEAR_INTEGRATION_AND_ENERGY_GPU_NO_SAMPLES(x, n_t, dt, t_0, orderTaylor, Adt, C_no_exp_all, C_no_exp_p, C_no_exp_k, C_no_exp_S, use_gpu::Bool, v_buf, vout_buf, term_buf)
    ## Allocate storage (CPU)
    E_all = zeros(n_t + 1)
    E_p = zeros(n_t + 1)
    E_k = zeros(n_t + 1)
    E_S = zeros(n_t + 1)
    t = zeros(n_t + 1)

    t[1] = t_0

    ## Copy initial vector into pre-allocated buffer
    copyto!(v_buf, x)
    copyto!(vout_buf, v_buf)

    ## Forward propagation via Taylor-series exponential (no allocation)
    for i in 1:n_t
        factorial_val = 1
        for j in 1:orderTaylor
            mul!(term_buf, Adt, v_buf, 1.0, 0.0)
            copyto!(v_buf, term_buf)
            factorial_val = j * factorial_val
            vout_buf .+= v_buf .* (1.0 / factorial_val)
        end
        copyto!(v_buf, vout_buf)

        # Store iteration data
        t[i+1] = t[i] + dt
        if use_gpu
            v_real = Array(real.(v_buf))
            v_real_dev = real.(v_buf)
        else
            v_real = real.(v_buf)
            v_real_dev = v_real
        end
        E_all[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_all * v_real_dev) : C_no_exp_all * v_real_dev))
        E_p[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_p * v_real_dev) : C_no_exp_p * v_real_dev))
        E_k[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_k * v_real_dev) : C_no_exp_k * v_real_dev))
        E_S[i+1] = real(dot(v_real, use_gpu ? Array(C_no_exp_S * v_real_dev) : C_no_exp_S * v_real_dev))
    end

    return (t, E_all, E_p, E_k, E_S)
end
