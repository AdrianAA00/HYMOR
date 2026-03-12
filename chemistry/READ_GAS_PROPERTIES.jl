using Printf
using Dates

"""
    READ_GAS_PROPERTIES(chemistry, filename)

Read gas properties from a data file into a Dict.

Parses a formatted gas properties data file and populates the chemistry Dict with
thermodynamic and species concentration data. Supports multiple file formats
(7, 16, 17, 18, or 19 columns) for backward compatibility.

# Arguments
- `chemistry::Dict{String,Any}`: Existing chemistry dict to augment
- `filename::String`: Path to the gas properties data file

# Returns
- `chemistry::Dict{String,Any}`: Updated dict containing gas properties

# Part of: Hypersonics Stability MATLAB Solver - Chemistry Module (Julia port)
"""
function READ_GAS_PROPERTIES(chemistry::Dict{String,Any}, filename::String)
    ## Validate input file
    if !isfile(filename)
        error("File \"$filename\" does not exist.")
    end

    try
        ## Parse header and data lines
        header_lines = String[]
        data_matrix  = Matrix{Float64}(undef, 0, 19)

        open(filename, "r") do fid
            for line in eachline(fid)
                # Skip empty lines
                stripped = strip(line)
                if isempty(stripped)
                    continue
                end

                if startswith(stripped, '#')
                    # Store header line
                    push!(header_lines, line)

                    # Extract specific metadata from header
                    if occursin("Planet:", line)
                        chemistry["planet"] = strip(split(line, "Planet:"; limit=2)[2])
                    elseif occursin("Gas model:", line)
                        chemistry["gas_model"] = strip(split(line, "Gas model:"; limit=2)[2])
                    elseif occursin("Number of valid data points:", line)
                        num_points_str = strip(split(line, "Number of valid data points:"; limit=2)[2])
                        chemistry["num_data_points"] = parse(Float64, num_points_str)
                    end
                else
                    # Parse numerical data line
                    values = parse.(Float64, split(stripped))
                    n_cols = length(values)

                    row = if n_cols == 19
                        # Full format: 9 properties + 10 species
                        reshape(values, 1, :)
                    elseif n_cols == 18
                        # Missing entropy: insert NaN at position 9
                        if size(data_matrix, 1) == 0
                            @warn "Format without entropy detected (18 columns). Entropy data will be NaN."
                        end
                        reshape(vcat(values[1:8], NaN, values[9:end]), 1, :)
                    elseif n_cols == 17
                        # Missing sound speed and entropy: insert NaN at positions 8-9
                        if size(data_matrix, 1) == 0
                            @warn "Legacy format detected (17 columns). Sound speed and entropy data will be NaN."
                        end
                        reshape(vcat(values[1:7], NaN, NaN, values[8:end]), 1, :)
                    elseif n_cols == 16
                        # Missing sound speed, entropy, and one species
                        reshape(vcat(values[1:7], NaN, NaN, values[8:end], NaN), 1, :)
                    elseif n_cols == 7
                        # Legacy format: basic properties only
                        if size(data_matrix, 1) == 0
                            @warn "Legacy format detected (7 columns). Sound speed, entropy, and molar data will be NaN."
                        end
                        reshape(vcat(values, fill(NaN, 12)), 1, :)
                    else
                        @warn "Skipping line with unexpected number of values ($n_cols): $line"
                        nothing
                    end

                    if row !== nothing
                        data_matrix = vcat(data_matrix, row)
                    end
                end
            end
        end

        chemistry["header"] = header_lines

        ## Validate parsed data
        if size(data_matrix, 1) == 0
            error("No valid numerical data found in file.")
        end

        ## Assign basic thermodynamic properties
        chemistry["rho"]        = data_matrix[:, 1]   # Density [kg/m^3]
        chemistry["e"]          = data_matrix[:, 2]   # Energy [J/kg]
        chemistry["T"]          = data_matrix[:, 3]   # Temperature [K]
        chemistry["gamma_star"] = data_matrix[:, 4]   # Heat capacity ratio [-]
        chemistry["cv_star"]    = data_matrix[:, 5]   # Specific heat [J/(kg*K)]
        chemistry["mu"]         = data_matrix[:, 6]   # Viscosity [Pa*s]
        chemistry["k"]          = data_matrix[:, 7]   # Thermal conductivity [W/(m*K)]
        chemistry["a"]          = data_matrix[:, 8]   # Sound speed [m/s]
        chemistry["s"]          = data_matrix[:, 9]   # Entropy [J/(kg*K)]

        ## Assign molar concentration data
        if size(data_matrix, 2) >= 19
            chemistry = assign_species_data(chemistry, data_matrix, 10, true)
        elseif size(data_matrix, 2) >= 18
            chemistry = assign_species_data(chemistry, data_matrix, 10, false)
        else
            chemistry["has_molar_data"] = false
            chemistry["species_list"] = String[]
        end

        ## Add metadata
        chemistry["num_points_actual"] = size(data_matrix, 1)
        chemistry["num_columns"]       = size(data_matrix, 2)
        chemistry["filename"]          = filename
        chemistry["date_read"]         = string(Dates.now())
        chemistry["has_sound_speed"]   = !all(isnan.(chemistry["a"]))
        chemistry["has_entropy"]       = !all(isnan.(chemistry["s"]))

        ## Display summary
        @printf("Successfully loaded gas properties data:\n")
        @printf("  File: %s\n", filename)
        if haskey(chemistry, "planet")
            @printf("  Planet: %s\n", chemistry["planet"])
        end
        if haskey(chemistry, "gas_model")
            @printf("  Gas Model: %s\n", chemistry["gas_model"])
        end
        @printf("  Data points read: %d\n", chemistry["num_points_actual"])
        @printf("  Columns: %d\n", chemistry["num_columns"])
        @printf("  Temperature range: %.1f - %.1f K\n", minimum(chemistry["T"]), maximum(chemistry["T"]))
        @printf("  Density range: %.2e - %.2e kg/m^3\n", minimum(chemistry["rho"]), maximum(chemistry["rho"]))

        if chemistry["has_sound_speed"]
            valid_a = chemistry["a"][.!isnan.(chemistry["a"])]
            @printf("  Sound speed range: %.1f - %.1f m/s\n", minimum(valid_a), maximum(valid_a))
        else
            @printf("  Sound speed data: Not available (legacy format)\n")
        end

        if chemistry["has_entropy"]
            valid_s = chemistry["s"][.!isnan.(chemistry["s"])]
            @printf("  Entropy range: %.1f - %.1f J/kg/K\n", minimum(valid_s), maximum(valid_s))
        else
            @printf("  Entropy data: Not available (legacy format)\n")
        end

        if chemistry["has_molar_data"]
            @printf("  Molar concentration data: Available\n")
            @printf("  Species included: %s\n", join(chemistry["species_list"], ", "))
            @printf("  Major species concentration ranges:\n")
            major_species = ["O2", "N2", "CO2", "H2"]
            for species in major_species
                if haskey(chemistry, species)
                    conc_data  = chemistry[species]
                    valid_data = conc_data[conc_data .> 0]
                    if !isempty(valid_data)
                        @printf("    %s: %.2e - %.2e mol\n", species, minimum(valid_data), maximum(valid_data))
                    end
                end
            end
        else
            @printf("  Molar concentration data: Not available (legacy format)\n")
        end

    catch ME
        rethrow(ME)
    end

    return chemistry
end

# ========================================================================
#  HELPER: Assign species concentration data to chemistry dict
# ========================================================================
"""
    assign_species_data(chemistry, data_matrix, col_start, has_N)

Extract species concentrations from data matrix.
"""
function assign_species_data(chemistry::Dict{String,Any}, data_matrix::Matrix{Float64}, col_start::Int, has_N::Bool)
    species_names = ["CO2", "H2", "O2", "N2", "CO", "NO", "C", "O", "H", "N"]

    # Assign log10 molar concentrations
    for i in 1:9
        field = "log10_" * species_names[i]
        chemistry[field] = data_matrix[:, col_start + i - 1]
    end
    if has_N
        chemistry["log10_N"] = data_matrix[:, col_start + 9]
    else
        chemistry["log10_N"] = fill(NaN, size(data_matrix, 1))
    end

    # Compute linear molar concentrations from log10 values
    for i in 1:length(species_names)
        log10_field  = "log10_" * species_names[i]
        linear_field = species_names[i]
        chemistry[linear_field] = 10.0 .^ chemistry[log10_field]
    end

    # Replace negligible concentrations (log10 <= -29.9) with zero
    active_species = has_N ? species_names : species_names[1:9]
    for species in active_species
        log10_field  = "log10_" * species
        linear_field = species
        low_conc_mask = chemistry[log10_field] .<= -29.9
        chemistry[linear_field][low_conc_mask] .= 0.0
    end

    # Store species metadata
    chemistry["species_list"]   = species_names
    chemistry["has_molar_data"] = true

    return chemistry
end
