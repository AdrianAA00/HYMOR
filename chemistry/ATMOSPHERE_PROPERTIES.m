function [rho, T] = ATMOSPHERE_PROPERTIES(planet, h)
% ATMOSPHERE_PROPERTIES  Compute atmospheric density and temperature at a given altitude.
%
%   [rho, T] = ATMOSPHERE_PROPERTIES(planet, h) returns the atmospheric
%   density and temperature for the specified planet and altitude using
%   standard atmosphere models.
%
% Inputs:
%   planet  - (string) Planet identifier: 'Earth' or 'Mars'
%   h       - (double) Altitude above the surface [m]
%
% Outputs:
%   rho     - (double) Atmospheric density [kg/m^3]
%   T       - (double) Atmospheric temperature [K]
%
% Models:
%   EARTH: U.S. Standard Atmosphere 1976 (NASA-TM-X-74335)
%          Five-layer model covering troposphere through stratopause (0-47 km).
%   MARS:  NASA Glenn Research Center Mars Atmosphere Model
%          Two-zone model with transition at 7 km altitude.
%          Reference: https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
%
% Examples:
%   [rho, T] = ATMOSPHERE_PROPERTIES('Earth', 10000)  % At 10 km altitude
%   [rho, T] = ATMOSPHERE_PROPERTIES('Mars', 5000)    % At 5 km altitude
%
% Part of: Hypersonics Stability MATLAB Solver - Chemistry Module

% Universal gas constant [J/(mol*K)]
R_UNIVERSAL = 8.31446;

% Select planet model
switch lower(planet)
    case 'earth'
        [rho, T] = EARTH_ATMOSPHERE(h, R_UNIVERSAL);
    case 'mars'
        [rho, T] = MARS_ATMOSPHERE(h);
    otherwise
        error('Planet must be ''Earth'' or ''Mars''');
end

end

%% ========================================================================
%  EARTH ATMOSPHERE MODEL (U.S. Standard Atmosphere 1976)
%  ========================================================================
function [rho, T] = EARTH_ATMOSPHERE(h, R_UNIVERSAL)
% EARTH_ATMOSPHERE  Compute density and temperature using the ISA model.

    % Earth constants
    g0 = 9.80665;           % Standard gravity [m/s^2]
    M  = 0.0289644;         % Molar mass of air [kg/mol]
    R  = R_UNIVERSAL / M;   % Specific gas constant [J/(kg*K)]

    % Layer boundaries and properties
    % [altitude (m), base temp (K), lapse rate (K/m)]
    layers = [
        0,     288.15, -0.0065;    % Troposphere
        11000, 216.65,  0.0;       % Tropopause / Lower Stratosphere
        20000, 216.65,  0.001;     % Stratosphere
        32000, 228.65,  0.0028;    % Upper Stratosphere
        47000, 270.65,  0.0;       % Stratopause
    ];

    % Base pressure at sea level [Pa]
    p0 = 101325;

    % Calculate temperature, pressure, and density
    T   = EARTH_TEMPERATURE(h, layers);
    p   = EARTH_PRESSURE(h, layers, p0, g0, R);
    rho = p * M / (R_UNIVERSAL * T);

end

%% ========================================================================
%  EARTH TEMPERATURE
%  ========================================================================
function T = EARTH_TEMPERATURE(h, layers)
% EARTH_TEMPERATURE  Calculate temperature at altitude h using the ISA model.

    n_layers = size(layers, 1);

    for i = 1:n_layers
        h_base = layers(i, 1);
        T_base = layers(i, 2);
        L      = layers(i, 3);

        if i < n_layers
            h_next = layers(i + 1, 1);
            if h <= h_next
                T = T_base + L * (h - h_base);
                return;
            end
        else
            % Last layer
            T = T_base + L * (h - h_base);
            return;
        end
    end

    % If beyond all layers, return last layer temperature
    T = layers(end, 2);

end

%% ========================================================================
%  EARTH PRESSURE
%  ========================================================================
function p = EARTH_PRESSURE(h, layers, p0, g0, R)
% EARTH_PRESSURE  Calculate pressure at altitude h using the ISA model.

    p = p0;
    h_prev   = 0;
    n_layers = size(layers, 1);

    for i = 1:n_layers
        h_base = layers(i, 1);
        T_base = layers(i, 2);
        L      = layers(i, 3);

        if i < n_layers
            h_next = layers(i + 1, 1);
        else
            h_next = h;
        end

        if h > h_base
            h_layer = min(h, h_next);
            dh      = h_layer - h_prev;

            if abs(L) > 1e-10  % Non-isothermal layer
                T_prev  = T_base + L * (h_prev - h_base);
                T_layer = T_base + L * (h_layer - h_base);
                p = p * (T_layer / T_prev)^(-g0 / (L * R));
            else  % Isothermal layer
                p = p * exp(-g0 * dh / (R * T_base));
            end

            h_prev = h_layer;

            if h <= h_next
                break;
            end
        end
    end

end

%% ========================================================================
%  MARS ATMOSPHERE MODEL (NASA Glenn Research Center)
%  ========================================================================
function [rho, T] = MARS_ATMOSPHERE(h)
% MARS_ATMOSPHERE  Compute density and temperature using the NASA GRC Mars model.
%   Two-zone model with temperature transition at 7 km altitude.

    % Mars constants
    h_transition = 7000;    % Transition altitude [m]
    R_mars       = 0.1921;  % Mars gas constant [kJ/(kg*K)]

    % Temperature calculation (two-zone model)
    if h <= h_transition
        T_celsius = -31.0 - 0.000998 * h;   % Lower atmosphere (0-7 km)
    else
        T_celsius = -23.4 - 0.00222 * h;    % Upper atmosphere (>7 km)
    end
    T = T_celsius + 273.15;  % Convert to Kelvin

    % Pressure calculation (exponential model) [kPa]
    p_kPa = 0.699 * exp(-0.00009 * h);

    % Density from equation of state: rho = p / (R * T)
    rho = p_kPa / (R_mars * T);  % kg/m^3

end
