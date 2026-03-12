function chemistry = READ_NON_EQUILIBRIUM_MODEL(chemistry, name)
% READ_NON_EQUILIBRIUM_MODEL - Reads non-equilibrium model coefficients from file
%
% Parses a structured text file containing non-equilibrium relaxation time
% model coefficients. Extracts metadata (planet, gas model, mechanism,
% temperature/density ranges) from comment lines and numeric coefficients
% for linear and quadratic fits of gamma* and cv* relaxation times.
%
% Syntax:
%   chemistry = READ_NON_EQUILIBRIUM_MODEL(chemistry, name)
%
% Inputs:
%   chemistry - Existing chemistry structure to augment with neq data
%   name      - Path to the non-equilibrium model file
%               (e.g., 'non_equi_model_Earth.txt')
%
% Outputs:
%   chemistry - Updated structure with chemistry.neq containing:
%               .planet          - Planet name string
%               .gas_model       - Gas model description
%               .mechanism       - Reaction mechanism
%               .T_range         - Temperature range string
%               .rho_range       - Density range string
%               .m_gamma         - Pressure exponent for gamma relaxation
%               .m_cv            - Pressure exponent for cv relaxation
%               .gamma.linear    - Linear fit coefficients (.a0, .a1, .R2)
%               .gamma.quadratic - Quadratic fit coefficients (.a0, .a1, .a2, .R2)
%               .cv.linear       - Linear fit coefficients (.a0, .a1, .R2)
%               .cv.quadratic    - Quadratic fit coefficients (.a0, .a1, .a2, .R2)
%               .model_equation  - Model equation string
%               .units           - Units structure (.tau, .P, .T)
%
% Notes:
%   The model equation is: ln(tau*P^m) = ln(T) + a0 + a1/T + a2/T^2
%   where tau is relaxation time (s), P is pressure (Pa), T is temperature (K).
%   The file must contain at least 16 numeric values in the expected order.
%
% Example:
%   chemistry = READ_NON_EQUILIBRIUM_MODEL(chemistry, 'non_equi_model_Earth.txt');
%   tau_gamma = exp(chemistry.neq.gamma.quadratic.a0 + ...
%                   chemistry.neq.gamma.quadratic.a1/T + ...
%                   chemistry.neq.gamma.quadratic.a2/T^2 + log(T)) / P^chemistry.neq.m_gamma;
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

    %% Initialize structure
    chemistry.neq = struct();

    %% Open and read file
    fid = fopen(name, 'r');
    if fid == -1
        error('Could not open file: %s', name);
    end

    %% Parse file contents
    lineNum = 0;
    valueIndex = 0;
    values = [];

    try
        % Read file line by line
        while ~feof(fid)
            line = fgetl(fid);
            lineNum = lineNum + 1;

            % Skip empty lines and comment lines starting with #
            if isempty(line) || startsWith(strtrim(line), '#')
                % Extract metadata from comments
                if contains(line, 'Planet:')
                    chemistry.neq.planet = strtrim(extractAfter(line, 'Planet:'));
                elseif contains(line, 'Gas Model:')
                    chemistry.neq.gas_model = strtrim(extractAfter(line, 'Gas Model:'));
                elseif contains(line, 'Mechanism:')
                    chemistry.neq.mechanism = strtrim(extractAfter(line, 'Mechanism:'));
                elseif contains(line, 'Temperature Range:')
                    chemistry.neq.T_range = strtrim(extractAfter(line, 'Temperature Range:'));
                elseif contains(line, 'Density Range:')
                    chemistry.neq.rho_range = strtrim(extractAfter(line, 'Density Range:'));
                end
                continue;
            end

            % Parse numeric values (format: value # comment)
            tokens = strsplit(line, '#');
            valueStr = strtrim(tokens{1});

            if ~isempty(valueStr)
                value = str2double(valueStr);
                if ~isnan(value)
                    valueIndex = valueIndex + 1;
                    values(valueIndex) = value;
                end
            end
        end

        fclose(fid);

        %% Assign parsed values to structure
        if length(values) >= 16
            % Pressure exponents
            chemistry.neq.m_gamma = values(1);
            chemistry.neq.m_cv = values(2);

            % Linear fit coefficients for gamma*
            chemistry.neq.gamma.linear.a0 = values(3);
            chemistry.neq.gamma.linear.a1 = values(4);
            chemistry.neq.gamma.linear.R2 = values(5);

            % Linear fit coefficients for cv*
            chemistry.neq.cv.linear.a0 = values(6);
            chemistry.neq.cv.linear.a1 = values(7);
            chemistry.neq.cv.linear.R2 = values(8);

            % Quadratic fit coefficients for gamma*
            chemistry.neq.gamma.quadratic.a0 = values(9);
            chemistry.neq.gamma.quadratic.a1 = values(10);
            chemistry.neq.gamma.quadratic.a2 = values(11);
            chemistry.neq.gamma.quadratic.R2 = values(12);

            % Quadratic fit coefficients for cv*
            chemistry.neq.cv.quadratic.a0 = values(13);
            chemistry.neq.cv.quadratic.a1 = values(14);
            chemistry.neq.cv.quadratic.a2 = values(15);
            chemistry.neq.cv.quadratic.R2 = values(16);

            % Store model equation and units
            chemistry.neq.model_equation = 'ln(tau*P^m) = ln(T) + a0 + a1/T + a2/T^2';
            chemistry.neq.units.tau = 's';
            chemistry.neq.units.P = 'Pa';
            chemistry.neq.units.T = 'K';
        else
            error('Insufficient data values found in file. Expected at least 16 values, found %d', length(values));
        end

    catch ME
        if fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end

    %% Display summary
    fprintf('Non-equilibrium model loaded successfully:\n');
    fprintf('  Planet: %s\n', chemistry.neq.planet);
    fprintf('  Mechanism: %s\n', chemistry.neq.mechanism);
    fprintf('  Temperature Range: %s\n', chemistry.neq.T_range);
    fprintf('  Pressure exponents: m_gamma = %.3f, m_cv = %.3f\n', ...
            chemistry.neq.m_gamma, chemistry.neq.m_cv);
    fprintf('  Linear fit R²: gamma = %.4f, cv = %.4f\n', ...
            chemistry.neq.gamma.linear.R2, chemistry.neq.cv.linear.R2);
    fprintf('  Quadratic fit R²: gamma = %.4f, cv = %.4f\n', ...
            chemistry.neq.gamma.quadratic.R2, chemistry.neq.cv.quadratic.R2);
end
