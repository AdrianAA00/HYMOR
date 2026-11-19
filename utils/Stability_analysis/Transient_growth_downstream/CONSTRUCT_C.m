function x = CONSTRUCT_C(x, prop, C, n_t, N_set_0)
% CONSTRUCT_C  Apply the transient growth operator C to a state vector.
%
%   x = CONSTRUCT_C(x, prop, C, n_t, N_set_0)
%
%   Computes C*x = exp(A'*t) * R'*M*R * exp(A*t) * x with the Taylor time
%   stepper TAYLOR_STEPS, for one vector or a block of them (the operator of
%   TRANSIENT_GROWTH_DOWNSTREAM and FREESTREAM_RECEPTIVITY).
%
%   Inputs:
%       x       - State vector or n x b block of them (CPU)
%       prop    - Time stepper from TAYLOR_PROPAGATOR (with the adjoint)
%       C       - Handle computing R'*M*R * x, from prop.matrix
%       n_t     - Number of time steps (scalar)
%       N_set_0 - Number of downstream DoFs zeroed at t = 0 and in the
%                 result (freestream receptivity; 0 otherwise)
%
%   Outputs:
%       x - Result of applying the operator C to the input (CPU)
%
% Part of: Hypersonics Stability MATLAB Solver - Transient Growth Downstream Module

    x(1:N_set_0, :) = 0;
    y = TAYLOR_STEPS(prop.load(x), prop.step, n_t);
    x = prop.unload(TAYLOR_STEPS(C(y), prop.step_adjoint, n_t));
    x(1:N_set_0, :) = 0;
end
