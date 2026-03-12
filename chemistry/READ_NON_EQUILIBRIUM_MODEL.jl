using Printf

"""
    READ_NON_EQUILIBRIUM_MODEL(chemistry, name)

Reads non-equilibrium model coefficients from file.

Parses a structured text file containing non-equilibrium relaxation time
model coefficients. Extracts metadata (planet, gas model, mechanism,
temperature/density ranges) from comment lines and numeric coefficients
for linear and quadratic fits of gamma* and cv* relaxation times.

# Arguments
- `chemistry::Dict{String,Any}`: Existing chemistry dict to augment with neq data
- `name::String`: Path to the non-equilibrium model file

# Returns
- `chemistry::Dict{String,Any}`: Updated dict with `chemistry["neq"]` containing model data

# Part of: Hypersonics Stability MATLAB Solver - Chemistry Module (Julia port)
"""
function READ_NON_EQUILIBRIUM_MODEL(chemistry::Dict{String,Any}, name::String)
    ## Initialize structure
    chemistry["neq"] = Dict{String,Any}()

    ## Open and read file
    if !isfile(name)
        error("Could not open file: $name")
    end

    ## Parse file contents
    lineNum = 0
    values = Float64[]

    try
        open(name, "r") do fid
            for line in eachline(fid)
                lineNum += 1

                # Skip empty lines and comment lines starting with #
                stripped = strip(line)
                if isempty(stripped) || startswith(stripped, '#')
                    # Extract metadata from comments
                    if occursin("Planet:", line)
                        chemistry["neq"]["planet"] = strip(split(line, "Planet:"; limit=2)[2])
                    elseif occursin("Gas Model:", line)
                        chemistry["neq"]["gas_model"] = strip(split(line, "Gas Model:"; limit=2)[2])
                    elseif occursin("Mechanism:", line)
                        chemistry["neq"]["mechanism"] = strip(split(line, "Mechanism:"; limit=2)[2])
                    elseif occursin("Temperature Range:", line)
                        chemistry["neq"]["T_range"] = strip(split(line, "Temperature Range:"; limit=2)[2])
                    elseif occursin("Density Range:", line)
                        chemistry["neq"]["rho_range"] = strip(split(line, "Density Range:"; limit=2)[2])
                    end
                    continue
                end

                # Parse numeric values (format: value # comment)
                tokens = split(line, '#')
                valueStr = strip(tokens[1])

                if !isempty(valueStr)
                    value = tryparse(Float64, valueStr)
                    if value !== nothing
                        push!(values, value)
                    end
                end
            end
        end

        ## Assign parsed values to structure
        if length(values) >= 16
            # Pressure exponents
            chemistry["neq"]["m_gamma"] = values[1]
            chemistry["neq"]["m_cv"] = values[2]

            # Linear fit coefficients for gamma*
            gamma_dict = Dict{String,Any}()
            gamma_linear = Dict{String,Any}()
            gamma_linear["a0"] = values[3]
            gamma_linear["a1"] = values[4]
            gamma_linear["R2"] = values[5]
            gamma_dict["linear"] = gamma_linear

            # Linear fit coefficients for cv*
            cv_dict = Dict{String,Any}()
            cv_linear = Dict{String,Any}()
            cv_linear["a0"] = values[6]
            cv_linear["a1"] = values[7]
            cv_linear["R2"] = values[8]
            cv_dict["linear"] = cv_linear

            # Quadratic fit coefficients for gamma*
            gamma_quadratic = Dict{String,Any}()
            gamma_quadratic["a0"] = values[9]
            gamma_quadratic["a1"] = values[10]
            gamma_quadratic["a2"] = values[11]
            gamma_quadratic["R2"] = values[12]
            gamma_dict["quadratic"] = gamma_quadratic

            # Quadratic fit coefficients for cv*
            cv_quadratic = Dict{String,Any}()
            cv_quadratic["a0"] = values[13]
            cv_quadratic["a1"] = values[14]
            cv_quadratic["a2"] = values[15]
            cv_quadratic["R2"] = values[16]
            cv_dict["quadratic"] = cv_quadratic

            chemistry["neq"]["gamma"] = gamma_dict
            chemistry["neq"]["cv"] = cv_dict

            # Store model equation and units
            chemistry["neq"]["model_equation"] = "ln(tau*P^m) = ln(T) + a0 + a1/T + a2/T^2"
            units = Dict{String,Any}()
            units["tau"] = "s"
            units["P"] = "Pa"
            units["T"] = "K"
            chemistry["neq"]["units"] = units
        else
            error("Insufficient data values found in file. Expected at least 16 values, found $(length(values))")
        end

    catch ME
        rethrow(ME)
    end

    ## Display summary
    @printf("Non-equilibrium model loaded successfully:\n")
    @printf("  Planet: %s\n", chemistry["neq"]["planet"])
    @printf("  Mechanism: %s\n", chemistry["neq"]["mechanism"])
    @printf("  Temperature Range: %s\n", chemistry["neq"]["T_range"])
    @printf("  Pressure exponents: m_gamma = %.3f, m_cv = %.3f\n",
            chemistry["neq"]["m_gamma"], chemistry["neq"]["m_cv"])
    @printf("  Linear fit R^2: gamma = %.4f, cv = %.4f\n",
            chemistry["neq"]["gamma"]["linear"]["R2"], chemistry["neq"]["cv"]["linear"]["R2"])
    @printf("  Quadratic fit R^2: gamma = %.4f, cv = %.4f\n",
            chemistry["neq"]["gamma"]["quadratic"]["R2"], chemistry["neq"]["cv"]["quadratic"]["R2"])

    return chemistry
end
