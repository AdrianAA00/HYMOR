# VERIFY_INPUT_DATA  Validate and augment the solver input data structure.
#
#   s = VERIFY_INPUT_DATA(s)
#
#   Performs comprehensive validation of all solver configuration fields,
#   checks types and allowed values, enforces consistency constraints
#   (e.g., periodic BCs must appear on both sides), and sets default values
#   for optional parameters that are not explicitly provided.
#
#   Inputs:
#       s - (Dict) Solver configuration structure.
#
#   Outputs:
#       s - (Dict) Validated structure with defaults set for optional
#           fields (numerical_dissipation, refinement, freestream
#           disturbance, stability_analysis.perturb_shock, etc.)
#
#   Notes:
#       - Chemistry-disabled mode sets defaults: chemistry_type = "None",
#         chemical_equilibrium = false, chemistry_composition = "None".
#       - Periodic BCs require matching on opposite boundaries and grid
#         dimensions divisible by 3.
#       - Supported boundary types: "circle", "lid_driven_cavity",
#         "channel", "MSL", "blunt_cone".
#       - Supported BC names: "inflow_subsonic", "inflow_supersonic",
#         "periodic", "shock", "outflow_supersonic", "outflow_subsonic",
#         "outflow_NRCBC", "no_slip_adiabatic", "no_slip_isothermal",
#         "symmetry".
#
# Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function VERIFY_INPUT_DATA(s::Dict{String,Any})

    # 1. Solver Options
    @assert haskey(s, "restart") "Missing field: restart"
    @assert isa(s["restart"], Bool) "restart must be a Bool"
    @assert haskey(s, "remesh") "Missing field: remesh"
    @assert isa(s["remesh"], Bool) "remesh must be a Bool"
    if s["restart"]
        @assert haskey(s, "restart_from_file") "Missing field: restart_from_file"
        @assert isa(s["restart_from_file"], Bool) "restart_from_file must be a Bool"
        if s["restart_from_file"]
            @assert haskey(s, "filename_restart") "Missing field: filename_restart"
            @assert isa(s["filename_restart"], String) "filename_restart must be a String"
        end
    end

    # 2. PDE dimension
    @assert haskey(s, "PDE_dimension") "Missing field: PDE_dimension"
    @assert s["PDE_dimension"] in ["2D", "3D-axisymmetric"] "PDE_dimension must be \"2D\" or \"3D-axisymmetric\""

    # 3. Chemistry
    @assert haskey(s, "chemistry") "Missing field: chemistry"
    @assert haskey(s["chemistry"], "is_chemistry_enabled") "Missing field: chemistry.is_chemistry_enabled"
    @assert isa(s["chemistry"]["is_chemistry_enabled"], Bool) "chemistry.is_chemistry_enabled must be a Bool"
    if s["chemistry"]["is_chemistry_enabled"]
        @assert haskey(s["chemistry"], "chemistry_type") "Missing field: chemistry.chemistry_type"
        @assert s["chemistry"]["chemistry_type"] in ["Frozen-RTV", "Frozen-RTVE", "Chemical-RTV", "Chemical-RTVE"] "Invalid chemistry_type"
        @assert haskey(s["chemistry"], "chemical_equilibrium") "Missing field: chemistry.chemical_equilibrium"
        @assert isa(s["chemistry"]["chemical_equilibrium"], Bool) "chemistry.chemical_equilibrium must be a Bool"
        @assert haskey(s["chemistry"], "non_equilibrium_model") "Missing field: chemistry.non_equilibrium_model"
        @assert s["chemistry"]["non_equilibrium_model"] in ["linear", "quadratic"] "Invalid non_equilibrium_model"
        @assert haskey(s["chemistry"], "chemistry_composition") "Missing field: chemistry.chemistry_composition"
        @assert s["chemistry"]["chemistry_composition"] in ["Earth", "Mars", "CO2"] "Invalid chemistry_composition"
    else
        s["chemistry"]["is_chemistry_enabled"]  = false
        s["chemistry"]["chemistry_type"]        = "None"
        s["chemistry"]["chemical_equilibrium"]  = false
        s["chemistry"]["chemistry_composition"] = "None"
        s["chemistry"]["planet"]                = "NONE"
        s["chemistry"]["gas_model"]             = "NONE"
    end

    # 4. Freestream
    @assert haskey(s, "freestream") "Missing field: freestream"
    @assert haskey(s["freestream"], "u") "Missing field: freestream.u"
    @assert isa(s["freestream"]["u"], Number) "freestream.u must be numeric"
    @assert haskey(s["freestream"], "v") "Missing field: freestream.v"
    @assert isa(s["freestream"]["v"], Number) "freestream.v must be numeric"
    if s["chemistry"]["is_chemistry_enabled"]
        @assert haskey(s["freestream"], "rho") "Missing field: freestream.rho"
        @assert isa(s["freestream"]["rho"], Number) && s["freestream"]["rho"] > 0 "freestream.rho must be positive numeric"
        @assert haskey(s["freestream"], "T") "Missing field: freestream.T"
        @assert isa(s["freestream"]["T"], Number) && s["freestream"]["T"] > 0 "freestream.T must be positive numeric"
    else
        @assert haskey(s["freestream"], "gamma") "Missing field: freestream.gamma"
        @assert isa(s["freestream"]["gamma"], Number) && s["freestream"]["gamma"] > 1 "freestream.gamma must be > 1"
        @assert haskey(s["freestream"], "Mach") "Missing field: freestream.Mach"
        @assert isa(s["freestream"]["Mach"], Number) && s["freestream"]["Mach"] > 0 "freestream.Mach must be positive numeric"
        @assert haskey(s["freestream"], "Pr") "Missing field: freestream.Pr"
        @assert isa(s["freestream"]["Pr"], Number) && s["freestream"]["Pr"] > 0 "freestream.Pr must be positive numeric"
    end

    if s["chemistry"]["is_chemistry_enabled"]
        has_L = haskey(s["freestream"], "L_ref")
        has_Re = haskey(s["freestream"], "Re")
        @assert xor(has_L, has_Re) "Exactly one of freestream.L_ref or freestream.Re must be specified, not both."

        if has_L
            @assert isa(s["freestream"]["L_ref"], Number) && s["freestream"]["L_ref"] > 0 "freestream.L_ref must be positive numeric"
        end
        if has_Re
            @assert isa(s["freestream"]["Re"], Number) && s["freestream"]["Re"] > 0 "freestream.Re must be positive numeric"
        end
    else
        @assert haskey(s["freestream"], "Re") "Missing field: freestream.Re"
        @assert isa(s["freestream"]["Re"], Number) && s["freestream"]["Re"] > 0 "freestream.Re must be positive numeric"
        @assert !haskey(s["freestream"], "L_ref") "freestream.L_ref should not be specified when chemistry is disabled. Specify freestream.Re instead."
        if !haskey(s["freestream"], "rho")
            s["freestream"]["rho"] = 1.0
        end
        if !haskey(s["freestream"], "T")
            s["freestream"]["T"] = 1.0
        end
    end

    # 4a. Coerce all numeric freestream parameters to Float64.
    #     Users may supply integer literals (e.g. T = 300, v = 0) in input
    #     files. Downstream code uses ::Float64 type assertions, so we
    #     convert here once, at the validation stage.
    for key in ("u", "v", "rho", "T", "Re", "L_ref",
                "gamma", "Mach", "Pr")
        if haskey(s["freestream"], key) && isa(s["freestream"][key], Number)
            s["freestream"][key] = Float64(s["freestream"][key])
        end
    end

    # 4b. Freestream Perturbation (optional - initialize to defaults if absent)
    if !haskey(s["freestream"], "disturbance") || isnothing(s["freestream"]["disturbance"]) || isempty(s["freestream"]["disturbance"])
        s["freestream"]["disturbance"] = Dict{String,Any}()
        s["freestream"]["disturbance"]["k_x"]       = 1.0
        s["freestream"]["disturbance"]["k_y"]       = 1.0
        s["freestream"]["disturbance"]["amplitude"] = [0.0, 0.0, 0.0, 0.0]
    else
        if !haskey(s["freestream"]["disturbance"], "k_x")
            s["freestream"]["disturbance"]["k_x"] = 1.0
        else
            @assert isa(s["freestream"]["disturbance"]["k_x"], Number) "freestream.disturbance.k_x must be numeric"
            s["freestream"]["disturbance"]["k_x"] = Float64(s["freestream"]["disturbance"]["k_x"])
        end
        if !haskey(s["freestream"]["disturbance"], "k_y")
            s["freestream"]["disturbance"]["k_y"] = 1.0
        else
            @assert isa(s["freestream"]["disturbance"]["k_y"], Number) "freestream.disturbance.k_y must be numeric"
            s["freestream"]["disturbance"]["k_y"] = Float64(s["freestream"]["disturbance"]["k_y"])
        end
        if !haskey(s["freestream"]["disturbance"], "amplitude")
            s["freestream"]["disturbance"]["amplitude"] = [0.0, 0.0, 0.0, 0.0]
        else
            @assert isa(s["freestream"]["disturbance"]["amplitude"], AbstractVector) && length(s["freestream"]["disturbance"]["amplitude"]) == 4 "freestream.disturbance.amplitude must be a numeric vector of length 4 [rho, u, v, rhoE]"
            s["freestream"]["disturbance"]["amplitude"] = Float64.(s["freestream"]["disturbance"]["amplitude"])
        end
    end

    # 5. Mesh
    @assert haskey(s, "mesh") "Missing field: mesh"
    @assert haskey(s["mesh"], "Nchi") "Missing field: mesh.Nchi"
    @assert isa(s["mesh"]["Nchi"], Number) && s["mesh"]["Nchi"] > 0 && mod(s["mesh"]["Nchi"], 1) == 0 "mesh.Nchi must be a positive integer"
    @assert haskey(s["mesh"], "Neta") "Missing field: mesh.Neta"
    @assert isa(s["mesh"]["Neta"], Number) && s["mesh"]["Neta"] > 0 && mod(s["mesh"]["Neta"], 1) == 0 "mesh.Neta must be a positive integer"
    s["mesh"]["Nchi"] = Int(s["mesh"]["Nchi"])
    s["mesh"]["Neta"] = Int(s["mesh"]["Neta"])

    # 6. Curvilinear Mapping
    @assert haskey(s, "curvilinear_mapping") "Missing field: curvilinear_mapping"
    @assert haskey(s["curvilinear_mapping"], "boundary_type") "Missing field: curvilinear_mapping.boundary_type"
    @assert s["curvilinear_mapping"]["boundary_type"] in ["circle", "lid_driven_cavity", "channel", "MSL", "blunt_cone"] "Invalid boundary_type"

    if !haskey(s["curvilinear_mapping"], "refinement_stagnation") || !haskey(s["curvilinear_mapping"]["refinement_stagnation"], "state")
        if !haskey(s["curvilinear_mapping"], "refinement_stagnation")
            s["curvilinear_mapping"]["refinement_stagnation"] = Dict{String,Any}()
        end
        s["curvilinear_mapping"]["refinement_stagnation"]["state"] = false
    end
    if !haskey(s["curvilinear_mapping"], "refinement_wall") || !haskey(s["curvilinear_mapping"]["refinement_wall"], "state")
        if !haskey(s["curvilinear_mapping"], "refinement_wall")
            s["curvilinear_mapping"]["refinement_wall"] = Dict{String,Any}()
        end
        s["curvilinear_mapping"]["refinement_wall"]["state"] = false
    end

    # 6b. Coerce curvilinear_mapping numeric parameters to Float64.
    for key in ("R", "dRe", "dRs", "circle_angle_extra",
                "channel_angle", "eta_refinement_power",
                "Lx", "Ly", "theta", "L")
        if haskey(s["curvilinear_mapping"], key) && isa(s["curvilinear_mapping"][key], Number)
            s["curvilinear_mapping"][key] = Float64(s["curvilinear_mapping"][key])
        end
    end
    # Refinement sub-dicts
    if haskey(s["curvilinear_mapping"], "refinement_stagnation")
        for key in ("BL_thickness", "intensity")
            if haskey(s["curvilinear_mapping"]["refinement_stagnation"], key) &&
               isa(s["curvilinear_mapping"]["refinement_stagnation"][key], Number)
                s["curvilinear_mapping"]["refinement_stagnation"][key] =
                    Float64(s["curvilinear_mapping"]["refinement_stagnation"][key])
            end
        end
    end

    # 7. Boundary Conditions
    @assert haskey(s, "boundary_conditions") "Missing field: boundary_conditions"
    valid_bcs = ["inflow_subsonic", "inflow_supersonic", "periodic", "shock", "outflow_supersonic", "outflow_subsonic", "outflow_NRCBC", "no_slip_adiabatic", "no_slip_isothermal", "symmetry"]
    @assert haskey(s["boundary_conditions"], "boundary_eta0") "Missing field: boundary_conditions.boundary_eta0"
    @assert s["boundary_conditions"]["boundary_eta0"]["name"] in valid_bcs "Invalid boundary_eta0"
    @assert haskey(s["boundary_conditions"], "boundary_eta1") "Missing field: boundary_conditions.boundary_eta1"
    @assert s["boundary_conditions"]["boundary_eta1"]["name"] in valid_bcs "Invalid boundary_eta1"
    @assert haskey(s["boundary_conditions"], "boundary_chi0") "Missing field: boundary_conditions.boundary_chi0"
    @assert s["boundary_conditions"]["boundary_chi0"]["name"] in valid_bcs "Invalid boundary_chi0"
    @assert haskey(s["boundary_conditions"], "boundary_chi1") "Missing field: boundary_conditions.boundary_chi1"
    @assert s["boundary_conditions"]["boundary_chi1"]["name"] in valid_bcs "Invalid boundary_chi1"

    # 7b. Periodic boundary condition consistency
    #     If one side of a direction is periodic, the opposite side must also be periodic.
    eta0_periodic = s["boundary_conditions"]["boundary_eta0"]["name"] == "periodic"
    eta1_periodic = s["boundary_conditions"]["boundary_eta1"]["name"] == "periodic"
    @assert eta0_periodic == eta1_periodic "Periodic boundary conditions must be set on both sides: boundary_eta0 and boundary_eta1 must both be periodic or both non-periodic."

    chi0_periodic = s["boundary_conditions"]["boundary_chi0"]["name"] == "periodic"
    chi1_periodic = s["boundary_conditions"]["boundary_chi1"]["name"] == "periodic"
    @assert chi0_periodic == chi1_periodic "Periodic boundary conditions must be set on both sides: boundary_chi0 and boundary_chi1 must both be periodic or both non-periodic."

    # 7c. When both eta boundaries are periodic, Neta must be a multiple of 3
    if eta0_periodic && eta1_periodic
        @assert mod(s["mesh"]["Neta"], 3) == 0 "When both boundary_eta0 and boundary_eta1 are periodic, mesh.Neta must be a multiple of 3."
    end
    if chi0_periodic && chi1_periodic
        @assert mod(s["mesh"]["Nchi"], 3) == 0 "When both boundary_chi0 and boundary_chi1 are periodic, mesh.Nchi must be a multiple of 3."
    end

    # 8. Time Integration
    @assert haskey(s, "time_integration") "Missing field: time_integration"
    @assert haskey(s["time_integration"], "N_iter") "Missing field: time_integration.N_iter"
    @assert isa(s["time_integration"]["N_iter"], Number) && s["time_integration"]["N_iter"] > 0 && mod(s["time_integration"]["N_iter"], 1) == 0 "time_integration.N_iter must be a positive integer"
    @assert haskey(s["time_integration"], "time_integrator") "Missing field: time_integration.time_integrator"
    @assert s["time_integration"]["time_integrator"] in ["Explicit_RK4", "Implicit_Euler"] "Invalid time_integrator"
    if s["time_integration"]["time_integrator"] == "Implicit_Euler"
        @assert haskey(s["time_integration"], "tolerance") "Missing field: time_integration.tolerance for Implicit_Euler"
        @assert isa(s["time_integration"]["tolerance"], Number) && s["time_integration"]["tolerance"] > 0 "time_integration.tolerance must be positive numeric"
        @assert haskey(s["time_integration"], "max_iter_implicit") "Missing field: time_integration.max_iter_implicit for Implicit_Euler"
        @assert isa(s["time_integration"]["max_iter_implicit"], Number) && s["time_integration"]["max_iter_implicit"] > 0 && mod(s["time_integration"]["max_iter_implicit"], 1) == 0 "time_integration.max_iter_implicit must be a positive integer"
        @assert haskey(s["time_integration"], "relax_factor") "Missing field: time_integration.relax_factor for Implicit_Euler"
        @assert isa(s["time_integration"]["relax_factor"], Number) && s["time_integration"]["relax_factor"] > 0 && s["time_integration"]["relax_factor"] <= 1 "time_integration.relax_factor must be between 0 and 1"
    end
    @assert haskey(s["time_integration"], "CFL") "Missing field: time_integration.CFL"
    @assert isa(s["time_integration"]["CFL"], Number) && s["time_integration"]["CFL"] > 0 "time_integration.CFL must be positive numeric"
    @assert haskey(s["time_integration"], "dt") "Missing field: time_integration.dt"
    @assert isa(s["time_integration"]["dt"], Number) && s["time_integration"]["dt"] > 0 "time_integration.dt must be positive numeric"
    @assert haskey(s["time_integration"], "max_dt") "Missing field: time_integration.max_dt"
    @assert isa(s["time_integration"]["max_dt"], Number) && s["time_integration"]["max_dt"] > 0 "time_integration.max_dt must be positive numeric"

    # 8b. Coerce time_integration numeric parameters to their proper types.
    s["time_integration"]["N_iter"] = Int(s["time_integration"]["N_iter"])
    for key in ("CFL", "dt", "max_dt")
        s["time_integration"][key] = Float64(s["time_integration"][key])
    end
    if s["time_integration"]["time_integrator"] == "Implicit_Euler"
        s["time_integration"]["max_iter_implicit"] = Int(s["time_integration"]["max_iter_implicit"])
        s["time_integration"]["tolerance"]         = Float64(s["time_integration"]["tolerance"])
        s["time_integration"]["relax_factor"]      = Float64(s["time_integration"]["relax_factor"])
    end

    # 9. Numerical Dissipation
    if !haskey(s, "numerical_dissipation")
        s["numerical_dissipation"] = Dict{String,Any}()
        s["numerical_dissipation"]["mu_rho"]   = 0.0
        s["numerical_dissipation"]["mu_rho_u"] = 0.0
        s["numerical_dissipation"]["mu_rho_v"] = 0.0
        s["numerical_dissipation"]["mu_rho_E"] = 0.0
        s["numerical_dissipation"]["rho_min"]  = 0.05
    else
        @assert haskey(s["numerical_dissipation"], "mu_rho") "Missing field: numerical_dissipation.mu_rho"
        @assert isa(s["numerical_dissipation"]["mu_rho"], Number) && s["numerical_dissipation"]["mu_rho"] >= 0 "numerical_dissipation.mu_rho must be non-negative numeric"
        @assert haskey(s["numerical_dissipation"], "mu_rho_u") "Missing field: numerical_dissipation.mu_rho_u"
        @assert isa(s["numerical_dissipation"]["mu_rho_u"], Number) && s["numerical_dissipation"]["mu_rho_u"] >= 0 "numerical_dissipation.mu_rho_u must be non-negative numeric"
        @assert haskey(s["numerical_dissipation"], "mu_rho_v") "Missing field: numerical_dissipation.mu_rho_v"
        @assert isa(s["numerical_dissipation"]["mu_rho_v"], Number) && s["numerical_dissipation"]["mu_rho_v"] >= 0 "numerical_dissipation.mu_rho_v must be non-negative numeric"
        @assert haskey(s["numerical_dissipation"], "mu_rho_E") "Missing field: numerical_dissipation.mu_rho_E"
        @assert isa(s["numerical_dissipation"]["mu_rho_E"], Number) && s["numerical_dissipation"]["mu_rho_E"] >= 0 "numerical_dissipation.mu_rho_E must be non-negative numeric"
        @assert haskey(s["numerical_dissipation"], "rho_min") "Missing field: numerical_dissipation.rho_min"
        @assert isa(s["numerical_dissipation"]["rho_min"], Number) && s["numerical_dissipation"]["rho_min"] >= 0 "numerical_dissipation.rho_min must be non-negative numeric"
        # Coerce to Float64
        for key in ("mu_rho", "mu_rho_u", "mu_rho_v", "mu_rho_E", "rho_min")
            s["numerical_dissipation"][key] = Float64(s["numerical_dissipation"][key])
        end
    end

    # 10. Shock Fitting
    @assert haskey(s, "shock") "Missing field: shock"
    @assert haskey(s["shock"], "enabled") "Missing field: shock.enabled"
    @assert isa(s["shock"]["enabled"], Bool) "shock.enabled must be a Bool"
    if s["shock"]["enabled"]
        @assert haskey(s["shock"], "feedback") "Missing field: shock.feedback"
        @assert isa(s["shock"]["feedback"], Bool) "shock.feedback must be a Bool"
        @assert haskey(s["shock"], "interpolate") "Missing field: shock.interpolate"
        @assert s["shock"]["interpolate"] in ["1st", "2nd", "3rd"] "Invalid shock.interpolate"
        @assert haskey(s["shock"], "initial_shock_dist") "Missing field: shock.initial_shock_dist"
        @assert isa(s["shock"]["initial_shock_dist"], Number) && s["shock"]["initial_shock_dist"] > 0 "shock.initial_shock_dist must be positive numeric"
        @assert haskey(s["shock"], "remesh_shock_distance") "Missing field: shock.remesh_shock_distance"
        @assert isa(s["shock"]["remesh_shock_distance"], Number) && s["shock"]["remesh_shock_distance"] > 0 "shock.remesh_shock_distance must be positive numeric"
        @assert haskey(s["shock"], "relaxation") "Missing field: shock.relaxation"
        @assert isa(s["shock"]["relaxation"], Number) && s["shock"]["relaxation"] >= 0 && s["shock"]["relaxation"] <= 1 "shock.relaxation must be between 0 and 1"
        @assert haskey(s["shock"], "formulation") "Missing field: shock.formulation"
        @assert s["shock"]["formulation"] in ["Lagrangian", "Eulerian"] "Invalid shock.formulation"
        @assert haskey(s["shock"], "fitting") "Missing field: shock.fitting"
        @assert s["shock"]["fitting"] in ["csaps"] "Invalid shock.fitting: only \"csaps\" is supported"
        @assert haskey(s["shock"], "spline_param") "Missing field: shock.spline_param"
        @assert isa(s["shock"]["spline_param"], Number) && s["shock"]["spline_param"] >= 0 && s["shock"]["spline_param"] <= 1 "shock.spline_param must be between 0 and 1"

        # 10b. Coerce shock numeric parameters to Float64.
        for key in ("initial_shock_dist", "remesh_shock_distance", "relaxation",
                     "spline_param", "initial_beta")
            if haskey(s["shock"], key) && isa(s["shock"][key], Number)
                s["shock"][key] = Float64(s["shock"][key])
            end
        end
    end

    # 11. Stability Analysis
    @assert haskey(s, "stability_analysis") "Missing field: stability_analysis"
    @assert haskey(s["stability_analysis"], "perturbation_magnitude") "Missing field: stability_analysis.perturbation_magnitude"
    @assert isa(s["stability_analysis"]["perturbation_magnitude"], Number) && s["stability_analysis"]["perturbation_magnitude"] > 0 "stability_analysis.perturbation_magnitude must be positive numeric"
    s["stability_analysis"]["perturbation_magnitude"] = Float64(s["stability_analysis"]["perturbation_magnitude"])
    @assert haskey(s["stability_analysis"], "eigenvalue_solver") "Missing field: stability_analysis.eigenvalue_solver"
    @assert s["stability_analysis"]["eigenvalue_solver"] in ["CPU_LU", "GPU_TIMESTEPPER_ARNOLDI"] "Invalid eigenvalue_solver"
    if !haskey(s["stability_analysis"], "perturb_shock")
        s["stability_analysis"]["perturb_shock"] = false
    else
        @assert isa(s["stability_analysis"]["perturb_shock"], Bool) "stability_analysis.perturb_shock must be a Bool"
    end

    # 12. Running Plot
    @assert haskey(s, "running_plot") "Missing field: running_plot"
    @assert haskey(s["running_plot"], "enabled") "Missing field: running_plot.enabled"
    @assert isa(s["running_plot"]["enabled"], Bool) "running_plot.enabled must be a Bool"
    if s["running_plot"]["enabled"]
        @assert haskey(s["running_plot"], "variable") "Missing field: running_plot.variable"
        @assert isa(s["running_plot"]["variable"], AbstractVector) "running_plot.variable must be a vector of strings"
        @assert haskey(s["running_plot"], "timesteps") "Missing field: running_plot.timesteps"
        @assert isa(s["running_plot"]["timesteps"], Number) && s["running_plot"]["timesteps"] > 0 && mod(s["running_plot"]["timesteps"], 1) == 0 "running_plot.timesteps must be a positive integer"
        s["running_plot"]["timesteps"] = Int(s["running_plot"]["timesteps"])
    end

    return s
end
