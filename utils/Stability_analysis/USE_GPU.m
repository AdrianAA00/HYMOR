function use_gpu = USE_GPU(s, A)
% USE_GPU  Choose and report the device of the linear time stepper.
%
%   use_gpu = USE_GPU(s, A)
%
%   s.stability_analysis.device is "gpu", "cpu" or "auto" (also when the
%   field is absent). "auto" takes the GPU when one is found and A has at
%   least 3e5 nonzeros. Below that a GPU Taylor step is bound by the
%   kernel-launch latency (about 0.25 ms per step) and the multithreaded
%   CPU product of MATVEC_HANDLE is faster (A100 vs 56 cores: equal near
%   2e5 nonzeros, GPU 2.5x faster at 1.7e6).
%
%   Inputs:
%       s - Solution structure (reads .stability_analysis.device if set).
%       A - Sparse system matrix the time stepper will apply.
%
%   Outputs:
%       use_gpu - true to run the time stepper on the GPU.
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis Module

    device = "auto";
    if isfield(s.stability_analysis, 'device')
        device = s.stability_analysis.device;
    end
    switch device
        case "gpu"
            use_gpu = true;
        case "cpu"
            use_gpu = false;
        otherwise
            use_gpu = gpuDeviceCount > 0 && nnz(A) >= 3e5;
    end

    if use_gpu
        disp(gpuDevice().Name + " selected.")
    else
        fprintf("Running on CPU (%d threads).\n", maxNumCompThreads)
    end
end
