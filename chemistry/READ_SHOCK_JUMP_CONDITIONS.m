function shock_jump_test = READ_SHOCK_JUMP_CONDITIONS(planet, chemistry_type)
% READ_SHOCK_JUMP_CONDITIONS - Reads shock jump condition files from Python solver
%
% Parses text files containing pre-computed shock jump conditions generated
% by the Python SD Toolbox solver. Extracts upstream freestream conditions
% and downstream post-shock state data (Mach, velocity, pressure, density,
% temperature) for validation of the MATLAB Rankine-Hugoniot solver.
%
% Syntax:
%   shock_jump_test = READ_SHOCK_JUMP_CONDITIONS(planet, chemistry_type)
%
% Inputs:
%   planet         - Planet name string: 'Earth' or 'Mars' (default: 'Earth')
%   chemistry_type - Chemistry model string (default: 'Chemical-RTVE'):
%                    'Chemical-RTV'  - Chemical equilibrium, R-T-V energy modes
%                    'Frozen-RTV'    - Frozen chemistry, R-T-V energy modes
%                    'Chemical-RTVE' - Chemical equilibrium, R-T-V-E energy modes
%                    'Frozen-RTVE'   - Frozen chemistry, R-T-V-E energy modes
%
% Outputs:
%   shock_jump_test - Structure containing:
%       .planet         - Planet name
%       .chemistry_type - Chemistry model used
%       .upstream       - Upstream conditions structure:
%           .P1    - Freestream pressure (Pa)
%           .rho1  - Freestream density (kg/m^3)
%           .T1    - Freestream temperature (K)
%           .e1    - Freestream specific internal energy (J/kg)
%           .a_s1  - Freestream speed of sound (m/s)
%       .data           - Post-shock data structure (field name derived from
%                         chemistry_type with '-' replaced by '_'):
%           .M1    - Upstream Mach numbers
%           .w1    - Upstream velocities (m/s)
%           .P2    - Post-shock pressures (Pa)
%           .rho2  - Post-shock densities (kg/m^3)
%           .T2    - Post-shock temperatures (K)
%
% Notes:
%   File path is constructed as:
%     <solver_dir>/chemistry/shock_jump_properties/shock_jump_conditions_<type>_<planet>.txt
%   The file header must contain a line with 'Upstream conditions:' followed
%   by key=value pairs for P1, rho1, T1, e1, and a_s1.
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

    %% Handle default arguments
    if nargin < 1
        planet = 'Earth';
    end

    if nargin < 2
        chemistry_type = 'Chemical-RTVE';
    end

    %% Validate chemistry type
    valid_types = {'Chemical-RTV', 'Frozen-RTV', 'Chemical-RTVE', 'Frozen-RTVE'};
    if ~ismember(chemistry_type, valid_types)
        error('Invalid chemistry_type. Must be one of: %s', strjoin(valid_types, ', '));
    end

    fprintf('Reading shock jump condition file for planet: %s, chemistry: %s\n', planet, chemistry_type);

    %% Initialize output structure
    shock_jump_test = struct();
    shock_jump_test.planet = planet;
    shock_jump_test.chemistry_type = chemistry_type;
    shock_jump_test.upstream = struct();
    shock_jump_test.data = struct();

    %% Construct filename and read data
    filename = sprintf(s.solver_dir + 'chemistry/shock_jump_properties/shock_jump_conditions_%s_%s.txt', chemistry_type, planet);

    fprintf('  Reading file: %s\n', filename);

    if exist(filename, 'file')
        try
            [data, upstream] = read_single_file(filename);

            % Store data with sanitized field name
            field_name = strrep(chemistry_type, '-', '_');
            shock_jump_test.data.(field_name) = data;

            % Store upstream conditions
            shock_jump_test.upstream = upstream;
            fprintf('    Upstream conditions: P1=%.2e Pa, rho1=%.3e kg/m^3, T1=%.1f K\n', ...
                    upstream.P1, upstream.rho1, upstream.T1);
            fprintf('                        e1=%.2e J/kg, a_s1=%.1f m/s\n', ...
                    upstream.e1, upstream.a_s1);

            fprintf('    Successfully read %d data points for %s\n', length(data.M1), chemistry_type);
        catch ME
            error('Failed to read file %s: %s', filename, ME.message);
        end
    else
        error('File not found: %s', filename);
    end

    fprintf('Successfully loaded shock jump conditions for %s model\n', chemistry_type);
end


%% ========================================================================
%  Local function: read_single_file
%  ========================================================================
function [data, upstream] = read_single_file(filename)
% READ_SINGLE_FILE - Parse a single shock jump condition file
%
% Reads the header to extract upstream conditions, then parses the
% numeric data block containing [M1, w1, P2, rho2, T2] columns.
%
% Inputs:
%   filename - Full path to the shock jump conditions text file
%
% Outputs:
%   data     - Structure with fields .M1, .w1, .P2, .rho2, .T2
%   upstream - Structure with fields .P1, .rho1, .T1, .e1, .a_s1

    data = struct();
    upstream = struct();

    %% Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    try
        %% Read header and extract upstream conditions
        upstream_found = false;
        while ~feof(fid)
            line = fgetl(fid);

            if ~ischar(line) || isempty(strtrim(line))
                continue;
            end

            % Check for upstream conditions line
            if contains(line, 'Upstream conditions:')
                upstream = parse_upstream_conditions(line);
                upstream_found = true;
                continue;
            end

            % First non-comment data line found
            if ~startsWith(strtrim(line), '#')
                break;
            end
        end

        if ~upstream_found
            warning('Upstream conditions not found in file header');
        end

        %% Read numerical data
        data_matrix = [];

        % Process current line if valid
        if ischar(line)
            values = str2num(line); %#ok<ST2NM>
            if length(values) == 5
                data_matrix = [data_matrix; values];
            end
        end

        % Continue reading rest of file
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line) && ~isempty(strtrim(line)) && ~startsWith(strtrim(line), '#')
                values = str2num(line); %#ok<ST2NM>
                if length(values) == 5
                    data_matrix = [data_matrix; values];
                end
            end
        end

        fclose(fid);

        %% Validate and assign data
        if isempty(data_matrix)
            error('No valid numerical data found in file');
        end

        data.M1   = data_matrix(:, 1);  % Mach number
        data.w1   = data_matrix(:, 2);  % Velocity (m/s)
        data.P2   = data_matrix(:, 3);  % Pressure (Pa)
        data.rho2 = data_matrix(:, 4);  % Density (kg/m^3)
        data.T2   = data_matrix(:, 5);  % Temperature (K)

    catch ME
        fclose(fid);
        rethrow(ME);
    end
end


%% ========================================================================
%  Local function: parse_upstream_conditions
%  ========================================================================
function upstream = parse_upstream_conditions(line)
% PARSE_UPSTREAM_CONDITIONS - Extract upstream state from header line
%
% Parses a line of the form:
%   # Upstream conditions: P1=1.013e+05 Pa, rho1=1.225e+00 kg/m^3, ...
%
% Inputs:
%   line - Header line string containing upstream conditions
%
% Outputs:
%   upstream - Structure with fields .P1, .rho1, .T1, .e1, .a_s1

    upstream = struct();

    % Extract P1
    P1_match = regexp(line, 'P1=([0-9.e+-]+)', 'tokens');
    if ~isempty(P1_match)
        upstream.P1 = str2double(P1_match{1}{1});
    end

    % Extract rho1
    rho1_match = regexp(line, 'rho1=([0-9.e+-]+)', 'tokens');
    if ~isempty(rho1_match)
        upstream.rho1 = str2double(rho1_match{1}{1});
    end

    % Extract T1
    T1_match = regexp(line, 'T1=([0-9.e+-]+)', 'tokens');
    if ~isempty(T1_match)
        upstream.T1 = str2double(T1_match{1}{1});
    end

    % Extract e1
    e1_match = regexp(line, 'e1=([0-9.e+-]+)', 'tokens');
    if ~isempty(e1_match)
        upstream.e1 = str2double(e1_match{1}{1});
    end

    % Extract a_s1
    as1_match = regexp(line, 'a_s1=([0-9.e+-]+)', 'tokens');
    if ~isempty(as1_match)
        upstream.a_s1 = str2double(as1_match{1}{1});
    end
end
