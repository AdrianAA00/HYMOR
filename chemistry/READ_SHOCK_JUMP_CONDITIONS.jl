using Printf

"""
    READ_SHOCK_JUMP_CONDITIONS(planet, chemistry_type)

Reads shock jump condition files from Python solver.

Parses text files containing pre-computed shock jump conditions generated
by the Python SD Toolbox solver. Extracts upstream freestream conditions
and downstream post-shock state data (Mach, velocity, pressure, density,
temperature) for validation of the MATLAB Rankine-Hugoniot solver.

# Arguments
- `planet::String`: Planet name string: "Earth" or "Mars" (default: "Earth")
- `chemistry_type::String`: Chemistry model string (default: "Chemical-RTVE")

# Returns
- `shock_jump_test::Dict{String,Any}`: Dict containing shock jump data

# Part of: Hypersonics Stability MATLAB Solver - Chemistry Module (Julia port)
"""
function READ_SHOCK_JUMP_CONDITIONS(planet::String="Earth", chemistry_type::String="Chemical-RTVE")
    ## Validate chemistry type
    valid_types = ["Chemical-RTV", "Frozen-RTV", "Chemical-RTVE", "Frozen-RTVE"]
    if !(chemistry_type in valid_types)
        error("Invalid chemistry_type. Must be one of: $(join(valid_types, ", "))")
    end

    @printf("Reading shock jump condition file for planet: %s, chemistry: %s\n", planet, chemistry_type)

    ## Initialize output structure
    shock_jump_test = Dict{String,Any}()
    shock_jump_test["planet"] = planet
    shock_jump_test["chemistry_type"] = chemistry_type
    shock_jump_test["upstream"] = Dict{String,Any}()
    shock_jump_test["data"] = Dict{String,Any}()

    ## Construct filename and read data
    filename = @sprintf("chemistry/shock_jump_properties/shock_jump_conditions_%s_%s.txt", chemistry_type, planet)

    @printf("  Reading file: %s\n", filename)

    if isfile(filename)
        try
            data, upstream = read_single_file(filename)

            # Store data with sanitized field name
            field_name = replace(chemistry_type, "-" => "_")
            shock_jump_test["data"][field_name] = data

            # Store upstream conditions
            shock_jump_test["upstream"] = upstream
            @printf("    Upstream conditions: P1=%.2e Pa, rho1=%.3e kg/m^3, T1=%.1f K\n",
                    upstream["P1"], upstream["rho1"], upstream["T1"])
            @printf("                        e1=%.2e J/kg, a_s1=%.1f m/s\n",
                    upstream["e1"], upstream["a_s1"])

            @printf("    Successfully read %d data points for %s\n", length(data["M1"]), chemistry_type)
        catch ME
            error("Failed to read file $filename: $(sprint(showerror, ME))")
        end
    else
        error("File not found: $filename")
    end

    @printf("Successfully loaded shock jump conditions for %s model\n", chemistry_type)
    return shock_jump_test
end


# ========================================================================
#  Local function: read_single_file
# ========================================================================
"""
    read_single_file(filename)

Parse a single shock jump condition file. Returns `(data, upstream)`.
"""
function read_single_file(filename::String)
    data = Dict{String,Any}()
    upstream = Dict{String,Any}()

    ## Open file
    if !isfile(filename)
        error("Could not open file: $filename")
    end

    try
        data_matrix = Matrix{Float64}(undef, 0, 5)
        upstream_found = false
        past_header = false

        open(filename, "r") do fid
            for line in eachline(fid)
                stripped = strip(line)
                if isempty(stripped)
                    continue
                end

                # Check for upstream conditions line
                if occursin("Upstream conditions:", line)
                    upstream = parse_upstream_conditions(line)
                    upstream_found = true
                    continue
                end

                # Check if it's a comment line
                if startswith(stripped, '#')
                    continue
                end

                # Parse numerical data line
                values = tryparse_line(stripped)
                if values !== nothing && length(values) == 5
                    data_matrix = vcat(data_matrix, reshape(values, 1, :))
                end
            end
        end

        if !upstream_found
            @warn "Upstream conditions not found in file header"
        end

        ## Validate and assign data
        if size(data_matrix, 1) == 0
            error("No valid numerical data found in file")
        end

        data["M1"]   = data_matrix[:, 1]  # Mach number
        data["w1"]   = data_matrix[:, 2]  # Velocity (m/s)
        data["P2"]   = data_matrix[:, 3]  # Pressure (Pa)
        data["rho2"] = data_matrix[:, 4]  # Density (kg/m^3)
        data["T2"]   = data_matrix[:, 5]  # Temperature (K)

    catch ME
        rethrow(ME)
    end

    return data, upstream
end

"""
    tryparse_line(line)

Try to parse a line of space-separated floats. Returns a Vector{Float64} or nothing.
"""
function tryparse_line(line::AbstractString)
    parts = split(strip(line))
    result = Float64[]
    for p in parts
        val = tryparse(Float64, p)
        if val === nothing
            return nothing
        end
        push!(result, val)
    end
    return result
end


# ========================================================================
#  Local function: parse_upstream_conditions
# ========================================================================
"""
    parse_upstream_conditions(line)

Extract upstream state from header line.
"""
function parse_upstream_conditions(line::String)
    upstream = Dict{String,Any}()

    # Extract P1
    P1_match = match(r"P1=([0-9.eE+-]+)", line)
    if P1_match !== nothing
        upstream["P1"] = parse(Float64, P1_match.captures[1])
    end

    # Extract rho1
    rho1_match = match(r"rho1=([0-9.eE+-]+)", line)
    if rho1_match !== nothing
        upstream["rho1"] = parse(Float64, rho1_match.captures[1])
    end

    # Extract T1
    T1_match = match(r"T1=([0-9.eE+-]+)", line)
    if T1_match !== nothing
        upstream["T1"] = parse(Float64, T1_match.captures[1])
    end

    # Extract e1
    e1_match = match(r"e1=([0-9.eE+-]+)", line)
    if e1_match !== nothing
        upstream["e1"] = parse(Float64, e1_match.captures[1])
    end

    # Extract a_s1
    as1_match = match(r"a_s1=([0-9.eE+-]+)", line)
    if as1_match !== nothing
        upstream["a_s1"] = parse(Float64, as1_match.captures[1])
    end

    return upstream
end
