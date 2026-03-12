function chemistry = READ_GAS_PROPERTIES(chemistry, filename)
% READ_GAS_PROPERTIES  Read gas properties from a data file into a structure.
%
%   chemistry = READ_GAS_PROPERTIES(chemistry, filename) parses a formatted
%   gas properties data file and populates the chemistry structure with
%   thermodynamic and species concentration data. Supports multiple file
%   formats (7, 16, 17, 18, or 19 columns) for backward compatibility.
%
% Inputs:
%   chemistry - (struct) Existing chemistry structure to augment
%   filename  - (string) Path to the gas properties data file
%
% Outputs:
%   chemistry - (struct) Updated structure containing gas properties:
%       .rho          - Density [kg/m^3]
%       .e            - Internal energy [J/kg]
%       .T            - Temperature [K]
%       .gamma_star   - Effective heat capacity ratio [-]
%       .cv_star      - Effective specific heat [J/(kg*K)]
%       .mu           - Dynamic viscosity [Pa*s]
%       .k            - Thermal conductivity [W/(m*K)]
%       .a            - Sound speed [m/s]
%       .s            - Entropy [J/(kg*K)]
%       .log10_<sp>   - Log10 molar concentrations for each species [mol]
%       .header       - Cell array of header comment lines
%       .planet       - Planet name (if specified in header)
%       .gas_model    - Gas model description (if specified in header)
%
% Notes:
%   - Header lines must begin with '#'.
%   - The file format is auto-detected based on the number of columns:
%       19 columns: full format (9 properties + 10 species)
%       18 columns: missing entropy (8 + 10 species)
%       17 columns: missing sound speed and entropy (7 + 10 species)
%       16 columns: missing sound speed, entropy, and one species (7 + 9)
%        7 columns: legacy format (basic properties only)
%   - Species list: CO2, H2, O2, N2, CO, NO, C, O, H, N
%   - Concentrations with log10 <= -29.9 are treated as effectively zero.
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

    %% Validate input file
    if ~exist(filename, 'file')
        error('File "%s" does not exist.', filename);
    end

    try
        fid = fopen(filename, 'r');
        if fid == -1
            error('Could not open file "%s".', filename);
        end

        %% Parse header and data lines
        header_lines = {};
        data_matrix  = [];

        while ~feof(fid)
            line = fgetl(fid);

            % Skip empty lines
            if ~ischar(line) || isempty(strtrim(line))
                continue;
            end

            if startsWith(strtrim(line), '#')
                % Store header line
                header_lines{end + 1} = line; %#ok<AGROW>

                % Extract specific metadata from header
                if contains(line, 'Planet:')
                    chemistry.planet = strtrim(extractAfter(line, 'Planet:'));
                elseif contains(line, 'Gas model:')
                    chemistry.gas_model = strtrim(extractAfter(line, 'Gas model:'));
                elseif contains(line, 'Number of valid data points:')
                    num_points_str = strtrim(extractAfter(line, 'Number of valid data points:'));
                    chemistry.num_data_points = str2double(num_points_str);
                end
            else
                % Parse numerical data line
                values = str2num(line); %#ok<ST2NM>
                n_cols = length(values);

                switch n_cols
                    case 19
                        % Full format: 9 properties + 10 species
                        data_matrix = [data_matrix; values]; %#ok<AGROW>
                    case 18
                        % Missing entropy: insert NaN at position 9
                        data_matrix = [data_matrix; values(1:8), NaN, values(9:end)]; %#ok<AGROW>
                        if size(data_matrix, 1) == 1
                            warning('Format without entropy detected (18 columns). Entropy data will be NaN.');
                        end
                    case 17
                        % Missing sound speed and entropy: insert NaN at positions 8-9
                        data_matrix = [data_matrix; values(1:7), NaN, NaN, values(8:end)]; %#ok<AGROW>
                        if size(data_matrix, 1) == 1
                            warning('Legacy format detected (17 columns). Sound speed and entropy data will be NaN.');
                        end
                    case 16
                        % Missing sound speed, entropy, and one species
                        data_matrix = [data_matrix; values(1:7), NaN, NaN, values(8:end), NaN]; %#ok<AGROW>
                    case 7
                        % Legacy format: basic properties only
                        data_matrix = [data_matrix; values, NaN(1, 12)]; %#ok<AGROW>
                        if size(data_matrix, 1) == 1
                            warning('Legacy format detected (7 columns). Sound speed, entropy, and molar data will be NaN.');
                        end
                    otherwise
                        warning('Skipping line with unexpected number of values (%d): %s', n_cols, line);
                end
            end
        end

        chemistry.header = header_lines;
        fclose(fid);

        %% Validate parsed data
        if isempty(data_matrix)
            error('No valid numerical data found in file.');
        end

        %% Assign basic thermodynamic properties
        chemistry.rho        = data_matrix(:, 1);   % Density [kg/m^3]
        chemistry.e          = data_matrix(:, 2);   % Energy [J/kg]
        chemistry.T          = data_matrix(:, 3);   % Temperature [K]
        chemistry.gamma_star = data_matrix(:, 4);   % Heat capacity ratio [-]
        chemistry.cv_star    = data_matrix(:, 5);   % Specific heat [J/(kg*K)]
        chemistry.mu         = data_matrix(:, 6);   % Viscosity [Pa*s]
        chemistry.k          = data_matrix(:, 7);   % Thermal conductivity [W/(m*K)]
        chemistry.a          = data_matrix(:, 8);   % Sound speed [m/s]
        chemistry.s          = data_matrix(:, 9);   % Entropy [J/(kg*K)]

        %% Assign molar concentration data
        if size(data_matrix, 2) >= 19
            chemistry = assign_species_data(chemistry, data_matrix, 10, true);
        elseif size(data_matrix, 2) >= 18
            chemistry = assign_species_data(chemistry, data_matrix, 10, false);
        else
            chemistry.has_molar_data = false;
            chemistry.species_list = {};
        end

        %% Add metadata
        chemistry.num_points_actual = size(data_matrix, 1);
        chemistry.num_columns       = size(data_matrix, 2);
        chemistry.filename          = filename;
        chemistry.date_read         = datestr(now);
        chemistry.has_sound_speed   = ~all(isnan(chemistry.a));
        chemistry.has_entropy       = ~all(isnan(chemistry.s));

        %% Display summary
        fprintf('Successfully loaded gas properties data:\n');
        fprintf('  File: %s\n', filename);
        if isfield(chemistry, 'planet')
            fprintf('  Planet: %s\n', chemistry.planet);
        end
        if isfield(chemistry, 'gas_model')
            fprintf('  Gas Model: %s\n', chemistry.gas_model);
        end
        fprintf('  Data points read: %d\n', chemistry.num_points_actual);
        fprintf('  Columns: %d\n', chemistry.num_columns);
        fprintf('  Temperature range: %.1f - %.1f K\n', min(chemistry.T), max(chemistry.T));
        fprintf('  Density range: %.2e - %.2e kg/m^3\n', min(chemistry.rho), max(chemistry.rho));

        if chemistry.has_sound_speed
            valid_a = chemistry.a(~isnan(chemistry.a));
            fprintf('  Sound speed range: %.1f - %.1f m/s\n', min(valid_a), max(valid_a));
        else
            fprintf('  Sound speed data: Not available (legacy format)\n');
        end

        if chemistry.has_entropy
            valid_s = chemistry.s(~isnan(chemistry.s));
            fprintf('  Entropy range: %.1f - %.1f J/kg/K\n', min(valid_s), max(valid_s));
        else
            fprintf('  Entropy data: Not available (legacy format)\n');
        end

        if chemistry.has_molar_data
            fprintf('  Molar concentration data: Available\n');
            fprintf('  Species included: %s\n', strjoin(chemistry.species_list, ', '));
            fprintf('  Major species concentration ranges:\n');
            major_species = {'O2', 'N2', 'CO2', 'H2'};
            for i = 1:length(major_species)
                species = major_species{i};
                if isfield(chemistry, species)
                    conc_data  = chemistry.(species);
                    valid_data = conc_data(conc_data > 0);
                    if ~isempty(valid_data)
                        fprintf('    %s: %.2e - %.2e mol\n', species, min(valid_data), max(valid_data));
                    end
                end
            end
        else
            fprintf('  Molar concentration data: Not available (legacy format)\n');
        end

    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end
end

%% ========================================================================
%  HELPER: Assign species concentration data to chemistry structure
%  ========================================================================
function chemistry = assign_species_data(chemistry, data_matrix, col_start, has_N)
% assign_species_data  Extract species concentrations from data matrix.

    species_names = {'CO2', 'H2', 'O2', 'N2', 'CO', 'NO', 'C', 'O', 'H', 'N'};

    % Assign log10 molar concentrations
    for i = 1:9
        field = ['log10_' species_names{i}];
        chemistry.(field) = data_matrix(:, col_start + i - 1);
    end
    if has_N
        chemistry.log10_N = data_matrix(:, col_start + 9);
    else
        chemistry.log10_N = NaN(size(data_matrix, 1), 1);
    end

    % Compute linear molar concentrations from log10 values
    for i = 1:length(species_names)
        log10_field  = ['log10_' species_names{i}];
        linear_field = species_names{i};
        chemistry.(linear_field) = 10.^chemistry.(log10_field);
    end

    % Replace negligible concentrations (log10 <= -29.9) with zero
    active_species = species_names;
    if ~has_N
        active_species = species_names(1:9);
    end
    for i = 1:length(active_species)
        log10_field  = ['log10_' active_species{i}];
        linear_field = active_species{i};
        low_conc_mask = chemistry.(log10_field) <= -29.9;
        chemistry.(linear_field)(low_conc_mask) = 0;
    end

    % Store species metadata
    chemistry.species_list   = species_names;
    chemistry.has_molar_data = true;
end
