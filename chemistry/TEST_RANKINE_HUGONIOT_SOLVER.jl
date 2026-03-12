using Printf
using Plots

"""
    TEST_RANKINE_HUGONIOT_SOLVER(chemistry, s)

Validate the Rankine-Hugoniot shock solver against reference Python / SD Toolbox data.

Runs three test cases:
1. Single moderate-shock solve (sanity check).
2. Mach-number sweep compared with SD Toolbox reference data.
3. Performance benchmark (100 simultaneous solves).

# Arguments
- `chemistry::Dict{String,Any}`: Fitted chemistry structure from FIT_CHEMISTRY
- `s::Dict{String,Any}`: Solution configuration structure

# Part of: Hypersonics Stability MATLAB Solver - Chemistry Module (Julia port)
"""
function TEST_RANKINE_HUGONIOT_SOLVER(chemistry::Dict{String,Any}, s::Dict{String,Any})
    ## Guard: chemistry must be enabled
    if !s["chemistry"]["is_chemistry_enabled"]
        @printf("No chemistry enabled...\n")
        s["chemistry"]["is_chemistry_enabled"]  = false
        s["chemistry"]["chemistry_type"]   = "None"
        s["chemistry"]["chemistry_composition"] = "None"
        return
    end

    if !haskey(chemistry, "eval_gamma_star")
        error("Chemistry structure must be fitted first. Run FIT_CHEMISTRY(chemistry)")
    end

    @printf("Testing vectorized non-linear solver...\n\n")

    ## Read reference shock-jump data
    @printf("=== Reading Python Shock Jump Conditions ===\n")
    planet = s["chemistry"]["chemistry_composition"]
    chemistry_type = s["chemistry"]["chemistry_type"]

    shock_jump_test = READ_SHOCK_JUMP_CONDITIONS(planet, chemistry_type)

    ## Test 1: single moderate shock
    @printf("\n=== Test Case 1: Single Moderate Shock ===\n")
    w_1 = 1000.0    # m/s
    rho_1 = 0.1     # kg/m^3
    e_1 = 5e5       # J/kg

    rho_2, e_2, w_2 = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(chemistry, [w_1], [rho_1], [e_1])

    @printf("Single case results: rho_2 = %.3e kg/m^3, e_2 = %.3e J/kg, w_2 = %.3e m/s\n\n",
        rho_2[1], e_2[1], w_2[1])

    ## Test 2: Mach-number sweep with Python data comparison
    @printf("=== Test Case 2: Mach Number Sweep with Python Data Comparison ===\n")

    if !isempty(shock_jump_test["upstream"])
        rho_1_plot = shock_jump_test["upstream"]["rho1"]
        e_1_plot   = shock_jump_test["upstream"]["e1"]
        a_s        = shock_jump_test["upstream"]["a_s1"]
        T_1_plot   = shock_jump_test["upstream"]["T1"]

        @printf("Using upstream conditions from Python files:\n")
        @printf("  rho_1 = %.3e kg/m^3\n", rho_1_plot)
        @printf("  e_1   = %.3e J/kg\n", e_1_plot)
        @printf("  T_1   = %.1f K\n", T_1_plot)
        @printf("  a_s   = %.1f m/s\n", a_s)
    else
        rho_1_plot = 0.1     # kg/m^3
        e_1_plot   = 2e5     # J/kg
        a_s        = 340.0   # m/s
        T_1_plot   = 300.0   # K

        @printf("Using default upstream conditions (Python data not available):\n")
        @printf("  rho_1 = %.3e kg/m^3\n", rho_1_plot)
        @printf("  e_1   = %.3e J/kg\n", e_1_plot)
        @printf("  T_1   = %.1f K\n", T_1_plot)
        @printf("  a_s   = %.1f m/s\n", a_s)
    end

    PLOT_TEST_RANKINE_HUGONIOT(chemistry, rho_1_plot, e_1_plot, a_s, T_1_plot, shock_jump_test)

    ## Test 3: performance benchmark
    @printf("\n=== Test Case 3: Performance Test ===\n")
    N_perf = 100
    w_1_perf = collect(range(500.0, 2000.0, length=N_perf))

    @printf("Performance test with %d points...\n", N_perf)
    elapsed_time = @elapsed begin
        rho_2_perf, e_2_perf, w_2_perf = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(
            chemistry, w_1_perf, rho_1_plot, e_1_plot)
    end

    @printf("Performance results:\n")
    @printf("  Total time: %.3f seconds\n", elapsed_time)
    @printf("  Time per solve: %.1f ms\n", 1000.0 * elapsed_time / N_perf)

    @printf("\n=== All Tests Completed Successfully ===\n")
end


# ========================================================================
#  Plotting helper
# ========================================================================
"""
    PLOT_TEST_RANKINE_HUGONIOT(chemistry, rho_1, e_1, a_s, T_1, shock_jump_test)

Plot density and temperature ratios vs Mach number with comparison
to Python / SD Toolbox reference data.
"""
function PLOT_TEST_RANKINE_HUGONIOT(chemistry::Dict{String,Any},
                                    rho_1::Float64, e_1_in::Float64,
                                    a_s::Float64, T_1::Float64,
                                    shock_jump_test::Dict{String,Any})
    ## Input validation
    if !haskey(chemistry, "eval_gamma_star")
        error("Chemistry structure must be fitted first. Run FIT_CHEMISTRY(chemistry)")
    end

    if !isfinite(T_1) || T_1 <= 0
        error(@sprintf("T_1 must be a positive finite number, got T_1 = %.2f", T_1))
    end

    if !isfinite(rho_1) || rho_1 <= 0
        error(@sprintf("rho_1 must be a positive finite number, got rho_1 = %.2e", rho_1))
    end

    # Recompute e_1 from the interpolation for consistency
    e_1 = chemistry["eval_e"](T_1, rho_1)

    ## Mach-number sweep
    M1_range = collect(range(1.5, 30.0, length=100))
    N_points = length(M1_range)
    w_1_range = M1_range .* a_s

    @printf("Computing MATLAB Rankine-Hugoniot relations for M1 = %.1f to %.1f...\n",
        minimum(M1_range), maximum(M1_range))
    @printf("Using %d Mach number points\n", N_points)
    @printf("State 1 conditions:\n")
    @printf("  rho_1 = %.3e kg/m^3\n", rho_1)
    @printf("  e_1   = %.3e J/kg\n", e_1)
    @printf("  T_1   = %.1f K\n", T_1)
    @printf("  a_s   = %.1f m/s\n", a_s)

    @printf("Solving vectorized Rankine-Hugoniot equations...\n")
    rho_2_range, e_2_range, w_2_range = SOLVE_RANKINE_HUGONIOT_CHEMISTRY(
        chemistry, w_1_range, rho_1, e_1)

    if any(isnan.(rho_2_range)) || any(isinf.(rho_2_range)) ||
       any(isnan.(e_2_range))   || any(isinf.(e_2_range))
        @warn "Solver produced NaN or Inf values. This may indicate convergence issues."
    end

    ## Compute ratios
    density_ratio = rho_2_range ./ rho_1

    eval_T = chemistry["eval_T"]
    T_2_range = eval_T(rho_2_range, e_2_range)
    if all(isnan.(T_2_range))
        error("All computed temperatures are NaN. Check chemistry interpolation data range.")
    end

    temperature_ratio = T_2_range ./ T_1

    # Filter invalid values for plotting
    valid_idx = isfinite.(density_ratio) .& isfinite.(temperature_ratio) .&
                (density_ratio .> 0) .& (temperature_ratio .> 0)

    if sum(valid_idx) == 0
        error("No valid solutions found. Check chemistry data and initial conditions.")
    end

    M1_plot = M1_range[valid_idx]
    density_ratio_plot = density_ratio[valid_idx]
    temperature_ratio_plot = temperature_ratio[valid_idx]

    @printf("Using %d valid points out of %d total points for plotting\n",
        sum(valid_idx), length(valid_idx))

    ## Create figure - Subplot 1: density ratio
    p1 = plot(M1_plot, density_ratio_plot, linewidth=2.5, color=:blue,
              label="Present Solver", xlabel="M_1", ylabel="rho_2/rho_1",
              xlims=(1.5, 35), grid=true, legend=:bottomright)

    if haskey(shock_jump_test, "chemistry_type") && haskey(shock_jump_test, "data")
        chemistry_type = shock_jump_test["chemistry_type"]
        field_name = replace(chemistry_type, "-" => "_")

        if haskey(shock_jump_test["data"], field_name)
            data = shock_jump_test["data"][field_name]
            @printf("Adding Python data to density plot for %s:\n", chemistry_type)

            python_density_ratio = data["rho2"] ./ rho_1
            idx_skip = 1:2:length(data["M1"])
            scatter!(p1, data["M1"][idx_skip], python_density_ratio[idx_skip],
                     markersize=3, color=:red, markerstrokecolor=:black,
                     markerstrokewidth=0.4, label="SD Toolbox $chemistry_type")

            @printf("  %s: %d points\n", chemistry_type, length(data["M1"]))
        end
    end

    density_max = maximum(density_ratio_plot)
    if isfinite(density_max) && density_max > 0
        ylims!(p1, (0, density_max * 1.1))
    else
        ylims!(p1, (0, 10))
    end

    ## Subplot 2: temperature ratio
    p2 = plot(M1_plot, temperature_ratio_plot, linewidth=2.5, color=:blue,
              label="Present Solver", xlabel="M_1", ylabel="T_2/T_1",
              xlims=(1.5, 35), grid=true, legend=:bottomright)

    if haskey(shock_jump_test, "chemistry_type") && haskey(shock_jump_test, "data")
        chemistry_type = shock_jump_test["chemistry_type"]
        field_name = replace(chemistry_type, "-" => "_")

        if haskey(shock_jump_test["data"], field_name)
            data = shock_jump_test["data"][field_name]
            @printf("Adding Python data to temperature plot for %s:\n", chemistry_type)

            python_temperature_ratio = data["T2"] ./ T_1
            idx_skip = 1:2:length(data["M1"])
            scatter!(p2, data["M1"][idx_skip], python_temperature_ratio[idx_skip],
                     markersize=3, color=:red, markerstrokecolor=:black,
                     markerstrokewidth=0.4, label="SD Toolbox $chemistry_type")

            @printf("  %s: %d points\n", chemistry_type, length(data["M1"]))
        end
    end

    temp_max = maximum(temperature_ratio_plot)
    if isfinite(temp_max) && temp_max > 0
        ylims!(p2, (0, temp_max * 1.1))
    else
        ylims!(p2, (0, 20))
    end

    ## Combine subplots
    combined_plot = plot(p1, p2, layout=(1, 2), size=(1000, 400))
    display(combined_plot)

    ## Print summary statistics
    @printf("\nMATLAB Results Summary:\n")
    @printf("Mach number range: %.1f - %.1f\n", minimum(M1_plot), maximum(M1_plot))
    @printf("Density ratio range: %.2f - %.2f\n", minimum(density_ratio_plot), maximum(density_ratio_plot))
    @printf("Temperature ratio range: %.2f - %.2f\n", minimum(temperature_ratio_plot), maximum(temperature_ratio_plot))
    @printf("Temperature range at state 2: %.1f - %.1f K\n", minimum(T_2_range), maximum(T_2_range))

    # Compare with Python reference data
    if haskey(shock_jump_test, "chemistry_type") && haskey(shock_jump_test, "data")
        chemistry_type = shock_jump_test["chemistry_type"]
        field_name = replace(chemistry_type, "-" => "_")

        if haskey(shock_jump_test["data"], field_name)
            data = shock_jump_test["data"][field_name]
            @printf("\nPython Data Summary for %s:\n", chemistry_type)
            @printf("  Mach range: %.2f - %.2f\n", minimum(data["M1"]), maximum(data["M1"]))
            @printf("  Density ratio range: %.2f - %.2f\n", minimum(data["rho2"]) / rho_1, maximum(data["rho2"]) / rho_1)
            @printf("  Temperature ratio range: %.2f - %.2f\n", minimum(data["T2"]) / T_1, maximum(data["T2"]) / T_1)
            @printf("  Temperature range: %.1f - %.1f K\n", minimum(data["T2"]), maximum(data["T2"]))
        end
    end
end
