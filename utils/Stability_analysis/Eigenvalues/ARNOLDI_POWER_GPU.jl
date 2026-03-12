using CUDA
using LinearAlgebra

"""
    ARNOLDI_POWER_GPU(x, Adt, timesteps, use_gpu, x_buf, vout_buf, term_buf, x_cpu_buf)

Matrix exponential action via Taylor expansion. Works on both GPU and CPU.

All buffers (x_buf, vout_buf, term_buf) are pre-allocated OUTSIDE the eigs
loop to avoid repeated memory allocation. Inside this function, only memory
copies are performed (no allocation).

# Arguments
- `x`: Input vector (CPU array or SubArray view from Arpack, length n).
- `Adt`: System matrix A*dt (CuArray/CuSparse or CPU sparse).
- `timesteps`: Number of time steps to advance (scalar).
- `use_gpu`: Boolean flag: true for GPU path, false for CPU path.
- `x_buf`: Pre-allocated buffer for input vector (CuArray or CPU Vector).
- `vout_buf`: Pre-allocated buffer for output accumulation.
- `term_buf`: Pre-allocated buffer for Taylor term.
- `x_cpu_buf`: Pre-allocated dense CPU buffer to collect Arpack SubArray views
  into contiguous memory before GPU transfer.

# Returns
- `x_out`: Result vector after applying the matrix exponential action (CPU array).

# Notes
- Uses a 5th-order Taylor expansion of exp(A*dt) per time step.
- Arpack frequently passes SubArray views rather than contiguous Vectors.
  copyto!(CuArray, SubArray) falls back to scalar indexing, so we first
  collect into x_cpu_buf (dense CPU), then transfer to GPU in one block copy.
- When use_gpu=false, all operations are performed on CPU.

Part of: Hypersonics Stability Julia Solver - Stability Analysis / Eigenvalues Module
"""
function ARNOLDI_POWER_GPU(x::AbstractVector, Adt, timesteps, use_gpu::Bool, x_buf, vout_buf, term_buf, x_cpu_buf)
    # Safely collect Arpack's view into a contiguous dense CPU array
    copyto!(x_cpu_buf, x)

    if use_gpu
        # Fast, block memory transfer from CPU dense to GPU dense
        copyto!(x_buf, x_cpu_buf)
        POWER_METHOD!(x_buf, Adt, timesteps, vout_buf, term_buf)
        x_out = Array(x_buf)
    else
        copyto!(x_buf, x_cpu_buf)
        POWER_METHOD!(x_buf, Adt, timesteps, vout_buf, term_buf)
        x_out = copy(x_buf)
    end

    return x_out
end

## ========================================================================
#  Taylor Expansion Time Stepper (uses pre-allocated buffers)
#  ========================================================================
"""
    POWER_METHOD!(v, A, timesteps, vout, term)

Advance a vector through multiple time steps using a truncated Taylor series
expansion of the matrix exponential. Uses pre-allocated buffers vout and term
to avoid memory allocation.
"""
function POWER_METHOD!(v, A, timesteps, vout, term)
    order = 5  # Taylor expansion order
    copyto!(vout, v)

    timesteps_val = timesteps isa AbstractArray ? Int(Array(timesteps)[]) : Int(timesteps)

    for i in 1:timesteps_val
        factorial_val = 1
        for j in 1:order
            mul!(term, A, v, 1.0, 0.0)
            copyto!(v, term)
            factorial_val = j * factorial_val
            vout .+= v .* (1.0 / factorial_val)
        end
        copyto!(v, vout)
    end
    copyto!(v, vout)
    return nothing
end
