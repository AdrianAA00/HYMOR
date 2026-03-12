function [s, solution_save, chemistry] = RESTART_FROM_FILE(s, solution_save, chemistry, Re, rho, T)
% RESTART_FROM_FILE  Load a saved simulation state from a MAT-file.
%
%   [s, solution_save, chemistry] = RESTART_FROM_FILE(s, solution_save, chemistry, Re, rho, T)
%
%   When restart-from-file mode is enabled, this function loads the solver
%   state from the MAT-file specified in s.filename_restart. After
%   loading, any freestream parameters (Reynolds number, density,
%   temperature) provided as optional arguments will override the values
%   stored in the file, enabling parameter sweeps from a common base state.
%
%   Inputs:
%       s      - (struct) Solver data structure with fields:
%                         .restart           - (logical) Global restart flag
%                         .restart_from_file - (logical) Load from file flag
%                         .filename_restart      - (char) Path to the MAT-file
%       solution_save - (struct) Previous s state (overwritten on load).
%       chemistry     - (struct) Chemistry model (overwritten on load).
%       Re            - (double, optional) Reynolds number override. Pass []
%                       to keep the value from the file.
%       rho           - (double, optional) Freestream density override [kg/m^3].
%                       Pass [] to keep the value from the file.
%       T             - (double, optional) Freestream temperature override [K].
%                       Pass [] to keep the value from the file.
%
%   Outputs:
%       s      - (struct) Loaded and optionally modified solver state.
%       solution_save - (struct) Loaded previous s state.
%       chemistry     - (struct) Loaded chemistry model.
%
%   Notes:
%       - If s.restart or s.restart_from_file is false,
%         the function returns immediately without modification.
%       - The MAT-file is expected to contain the variables: s,
%         solution_save, and chemistry.
%
%   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Handle optional arguments
    if nargin < 6
        T = [];
    end
    if nargin < 5
        rho = [];
    end
    if nargin < 4
        Re = [];
    end

    %% Load state from file if restart-from-file is active
    if s.restart && s.restart_from_file
        load(s.filename_restart)
        disp("Restarting from file: " + s.filename_restart)
    else
        return
    end

    %% Override freestream parameters if provided
    disp(" ")
    disp("-----------------------")
    disp("Forced input parameters")

    if ~isempty(Re)
        s.freestream.Re = Re;
        disp("Re = " + s.freestream.Re)
    end
    if ~isempty(rho)
        s.freestream.rho = rho;
        disp("rho = " + s.freestream.rho)
    end
    if ~isempty(T)
        s.freestream.T = T;
        disp("T = " + s.freestream.T)
    end

    disp("-----------------------")
    disp(" ")
end
