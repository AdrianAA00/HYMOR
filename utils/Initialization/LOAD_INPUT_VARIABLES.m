function LOAD_INPUT_VARIABLES(filename)
% LOAD_INPUT_VARIABLES  Execute an input parameter script in the caller workspace.
%
%   LOAD_INPUT_VARIABLES(filename)
%
%   Reads the contents of the specified MATLAB script file and executes it
%   in the caller's workspace using evalin. This allows input parameter
%   files to define variables directly in the workspace of the calling
%   function or script without requiring explicit output assignments.
%
%   Inputs:
%       filename - (char) Path to the MATLAB script file (.m) containing
%                  variable definitions and configuration parameters.
%
%   Outputs:
%       (none)   - Variables are created/modified in the caller's workspace.
%
%   Notes:
%       - The input file must be a valid MATLAB script (not a function).
%       - An error is thrown if the specified file does not exist.
%       - Use with caution: evalin can overwrite existing variables in the
%         caller's workspace without warning.
%
%   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

    %% Validate that the file exists
    if ~exist(filename, 'file')
        error('File does not exist: %s', filename);
    end

    %% Read and execute the file contents in the caller workspace
    file_content = fileread(filename);
    evalin('caller', file_content);

end
