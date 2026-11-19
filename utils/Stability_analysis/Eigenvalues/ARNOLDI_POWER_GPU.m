function x = ARNOLDI_POWER_GPU(x, prop, timesteps)
% ARNOLDI_POWER_GPU  Action of the Taylor time stepper exp(A*dt*timesteps).
%
%   x = ARNOLDI_POWER_GPU(x, prop, timesteps) is the eigs callback of
%   EIGENVALUES (GPU_TIMESTEPPER_ARNOLDI): it moves x to the device,
%   advances it with TAYLOR_STEPS and returns it on the host.
%
%   Inputs:
%       x         - Input vector (CPU array, length n).
%       prop      - Time stepper from TAYLOR_PROPAGATOR.
%       timesteps - Number of time steps to advance (scalar).
%
%   Outputs:
%       x - Result vector after applying the matrix exponential action
%           (CPU array, length n).
%
% Part of: Hypersonics Stability MATLAB Solver - Stability Analysis / Eigenvalues Module

    x = prop.unload(TAYLOR_STEPS(prop.load(x), prop.step, timesteps));
end
