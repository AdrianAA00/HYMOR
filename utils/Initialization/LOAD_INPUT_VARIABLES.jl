# LOAD_INPUT_VARIABLES  Execute an input parameter script in the caller workspace.
#
#   LOAD_INPUT_VARIABLES(filename)
#
#   Reads the contents of the specified Julia script file and executes it
#   using include. This allows input parameter files to define variables
#   directly.
#
#   Inputs:
#       filename - (String) Path to the Julia script file (.jl) containing
#                  variable definitions and configuration parameters.
#
#   Outputs:
#       (none)   - Variables are created/modified via include.
#
#   Notes:
#       - The input file must be a valid Julia script.
#       - An error is thrown if the specified file does not exist.
#
#   Part of: Hypersonics Stability MATLAB Solver - Initialization Module

function LOAD_INPUT_VARIABLES(filename::String)

    ## Validate that the file exists
    if !isfile(filename)
        error("File does not exist: $filename")
    end

    ## Read and execute the file contents in Main scope so that variables
    ## (e.g. solution, solution_save) are visible to the calling script.
    Base.include(Main, filename)

end
