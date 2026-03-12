using CUDA
using LinearAlgebra

"""
    CONSTRUCT_C(x, A_dt, At_dt, C_no_exp, n_t, freestream_disturbances, N_set_0, use_gpu, x_buf, vout_buf, term_buf, x_cpu_buf)

Apply the transient growth operator C to a state vector. Works on both GPU and CPU.

All buffers (x_buf, vout_buf, term_buf) are pre-allocated OUTSIDE the eigs
loop to avoid repeated memory allocation. Inside this function, only memory
copies are performed (no allocation).

# Arguments
- `x`: State vector (column vector, may be a SubArray view from Arpack).
- `A_dt`: System Jacobian matrix A*dt (CuArray/CuSparse or CPU sparse).
- `At_dt`: Transpose (A*dt)' (CuArray/CuSparse or CPU sparse).
- `C_no_exp`: Energy weight matrix R'*M*R (CuArray/CuSparse or CPU sparse).
- `n_t`: Number of time steps (scalar).
- `freestream_disturbances`: Flag: true to zero out downstream DoFs at t=0.
- `N_set_0`: Number of downstream DoFs to zero when freestream mode is active.
- `use_gpu`: Boolean flag: true for GPU path, false for CPU path.
- `x_buf`: Pre-allocated buffer for input vector.
- `vout_buf`: Pre-allocated buffer for output accumulation.
- `term_buf`: Pre-allocated buffer for Taylor term.
- `x_cpu_buf`: Pre-allocated dense CPU buffer to collect Arpack SubArray views
  into contiguous memory before GPU transfer.

# Returns
- `x`: Result of applying the operator C to the input vector.

Part of: Hypersonics Stability Julia Solver - Transient Growth Downstream Module
"""
function CONSTRUCT_C(x, A_dt, At_dt, C_no_exp, n_t, freestream_disturbances::Bool, N_set_0::Int, use_gpu::Bool, x_buf, vout_buf, term_buf, x_cpu_buf)
    # Safely collect Arpack's view into a contiguous dense CPU array
    copyto!(x_cpu_buf, x)

    ## Zero downstream disturbances if freestream mode
    if freestream_disturbances
        x_cpu_buf[1:N_set_0] .= 0
    end

    if use_gpu
        ## GPU path: fast block copy from contiguous CPU to GPU
        copyto!(x_buf, x_cpu_buf)
        EXPONENTIAL!(x_buf, A_dt, At_dt, C_no_exp, n_t, vout_buf, term_buf)
        x = Array(x_buf)
    else
        ## CPU path: copy into pre-allocated buffer (no allocation)
        copyto!(x_buf, x_cpu_buf)
        EXPONENTIAL!(x_buf, A_dt, At_dt, C_no_exp, n_t, vout_buf, term_buf)
        x = copy(x_buf)
    end

    return x
end


"""
    EXPONENTIAL!(v, A, At, C_no_exp, timesteps, vout, term)

Compute exp(A'*t) * C_no_exp * exp(A*t) * v via Taylor series.
Uses pre-allocated buffers vout and term to avoid memory allocation.
Result is stored in v.

# Arguments
- `v`: Input state vector (modified in-place with result).
- `A`: System matrix (A*dt already incorporated).
- `At`: Transpose of system matrix (A*dt)'.
- `C_no_exp`: Energy weight matrix (R'*M*R).
- `timesteps`: Number of sub-steps for the exponential.
- `vout`: Pre-allocated buffer for output accumulation.
- `term`: Pre-allocated buffer for Taylor term.
"""
function EXPONENTIAL!(v, A, At, C_no_exp, timesteps, vout, term)
    order = 5  # Fixed order for Taylor expansion

    timesteps_val = timesteps isa AbstractArray ? Int(Array(timesteps)[]) : Int(timesteps)

    copyto!(vout, v)
    for i in 1:timesteps_val
        copyto!(term, vout)
        for j in 1:order
            mul!(v, A, term, 1.0, 0.0)     # v = A * term  (reuse v as temp)
            copyto!(term, v)             # term = v
            term .*= (1.0 / j)           # term = (A*dt)^j / j!
            vout .+= term
        end
    end

    ## Apply energy weight: v = R'*M*R * exp(A*t) * v
    mul!(v, C_no_exp, vout, 1.0, 0.0)
    copyto!(vout, v)

    ## Adjoint propagation: vout = exp(A'*t) * R'*M*R * exp(A*t) * v
    for i in 1:timesteps_val
        copyto!(term, vout)
        for j in 1:order
            mul!(v, At, term, 1.0, 0.0)    # v = At * term  (reuse v as temp)
            copyto!(term, v)     # term = v
            term .*= (1.0 / j)           # term = (At*dt)^j / j!
            vout .+= term
        end
    end

    copyto!(v, vout)
    return nothing
end
