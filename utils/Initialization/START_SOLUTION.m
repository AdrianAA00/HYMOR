function [s] = START_SOLUTION(s, chemistry)
% START_SOLUTION  Build and initialize a fresh s from scratch.
%   Allocates field arrays, generates the boundary and mesh, sets the
%   initial flow field, and (when a shock is present) initializes the
%   shock line, post-shock state, and shock boundary conditions.
%   Finishes by extending the upstream flow and applying all boundary
%   conditions.
%
%   [s] = START_SOLUTION(s, chemistry)
%
%   Inputs:
%       s  - Solution struct with configuration parameters
%                   (Nx, Ny, boundary_type, shock flag, etc.)
%       chemistry - Chemistry model struct for thermodynamic evaluations
%
%   Outputs:
%       s  - Fully initialized s struct ready for time
%                   integration, including mesh coordinates, conserved
%                   variables, thermodynamic fields, and boundary data
%
%   Notes:
%       - Time and iteration counter are reset to zero.
%       - Mesh is generated without disturbances for the initial setup.
%       - The shock initialization branch calls a sequence of sub-routines:
%         SHOCK_LINE_INITIALIZATION, LEAST_SQUARES_SHOCK_POINTS,
%         DETECT_CELLS_SHOCKED, COMPUTE_BETA, and
%         INITIAL_SOLUTION_POST_SHOCK.
%       - Chemistry equilibrium and thermodynamic properties are updated
%         after both the pre-shock and post-shock initialization stages.
%
% Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Reset time and iteration counter
    s.time_integration.t = 0;
    s.time_integration.iter = 0;
    s.linearize = false;

    %% Allocate arrays and generate mesh
    s = VARIABLES_INITIALIZATION(s);
    disturbances = false;
    s = GENERATE_MESH(s, disturbances);

    %% Set initial flow field and update thermodynamic state
    s = INITIAL_SOLUTION(s, chemistry);
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);

    %% Shock initialization (if shock is present)
    if s.shock.enabled
        s = SHOCK_LINE_INITIALIZATION(s);
        s = LEAST_SQUARES_SHOCK_POINTS(s);
        s = DETECT_CELLS_SHOCKED(s);
        s = COMPUTE_BETA(s);
        s = INITIAL_SOLUTION_POST_SHOCK(s, chemistry);
        % s = INITIAL_SOLUTION_POST_SHOCK_MODIFIED(s, chemistry);
        s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
        s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
        s = UPDATE_SOUND_SPEED(s, chemistry);
        s = UPDATE_SHOCK_BC(s, chemistry);
        s = UPDATE_FIELD_UPSTREAM(s);
    end

    %% Extend upstream flow and apply boundary conditions
    s = UPDATE_CHEMISTRY_EQUILIBRIUM(s, chemistry);
    s = UPDATE_THERMODYNAMIC_PROPERTIES(s, chemistry);
    s = APPLY_BOUNDARY_CONDITIONS(s, chemistry);
    s = UPDATE_SOUND_SPEED(s, chemistry);
end
