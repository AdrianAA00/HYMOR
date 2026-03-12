using LinearAlgebra
using SparseArrays
using Arpack
using Printf

"""
    FREESTREAM_RECEPTIVITY(A_, s, n_modes, T, w_infty)

Compute optimal transient growth modes including freestream disturbances.

Solves the generalized eigenvalue problem for transient energy growth
in the coupled downstream-freestream system. Uses a matrix-free
approach with eigs() to find the n_modes largest-gain eigenmodes at
the specified time horizon T. Automatically detects GPU availability.

# Arguments
- `A_`: Extended system matrix (sparse, from LINEARIZE_L_).
- `s`: Solution Dict{String,Any} with flow fields, mesh data, and configuration flags.
- `n_modes`: Number of leading transient growth modes to compute.
- `T`: Target time horizon for the optimization.
- `w_infty`: Vector of freestream disturbance frequencies.

# Returns
- `V`: Matrix of optimal perturbation eigenvectors (columns).
- `D`: Diagonal matrix of corresponding energy gains (eigenvalues).
- `T_opt`: Actual time horizon used (adjusted to integer time steps).

# Notes
- All memory is pre-allocated outside eigs. Only copies inside callback.
- GPU is synchronized and memory verified deallocated before returning.
- Falls back to CPU if no CUDA device is found.
- Downstream DoFs are zeroed at t=0 so that only freestream forcing
  drives the response (freestream_disturbances = true).

Part of: Hypersonics Stability Julia Solver - Transient Growth Freestream Module
"""
function FREESTREAM_RECEPTIVITY(A_, s::Dict{String,Any}, n_modes::Int, T, w_infty)
    k = n_modes        # Maximum of modes to compute
    freestream_disturbances = true

    ## Compute time-stepping parameters
    s = CFL_TIMESTEP(s)
    t_temp = s["time_integration"]["dt"]
    n_t = round(Int, T / t_temp)
    dt = T / n_t
    T_opt = n_t * dt

    ## Construct operator D_ and C_ (without exponentials)
    R_ = CONSTRUCT_R_(s, w_infty)
    norms_1 = [1, 1, 1, 1]
    norms_2 = [1, 1, 1, 1]
    M_ = CONSTRUCT_M_(s, norms_1, w_infty)
    scaling_non_temporal = true  # Fixed now
    M_infty_ = CONSTRUCT_M_INFTY_(s, norms_2, T, w_infty, scaling_non_temporal)
    C_no_exp_ = R_' * M_ * R_
    D_ = R_' * M_infty_ * R_

    ## Eigenvalue computation
    println()
    println("-----------------------------------")
    println("      Freestream receptivity       ")
    println("-----------------------------------")
    println()
    @printf("Number timesteps = %f \n", n_t)

    ## Detect GPU availability
    use_gpu = CUDA.functional()
    N_set_0 = 4 * s["mesh"]["Nchi"] * s["mesh"]["Neta"]
    n = size(A_, 1)

    if use_gpu
        gpu = CUDA.device()
        println("$(CUDA.name(gpu)) selected.")
    else
        println("No CUDA device found. Running on CPU.")
    end

    ## Pre-allocate ALL memory outside eigs (matrices + work buffers)
    A_dt = A_ * dt
    At_dt = sparse(A_dt')

    if use_gpu
        gA_dt      = CuSparseMatrixCSC(A_dt)
        gAt_dt     = CuSparseMatrixCSC(At_dt)
        gC_no_exp_ = CuSparseMatrixCSC(C_no_exp_)
        x_buf      = CUDA.zeros(ComplexF64, n)
        vout_buf   = CUDA.zeros(ComplexF64, n)
        term_buf   = CUDA.zeros(ComplexF64, n)
        x_cpu_buf  = zeros(ComplexF64, n)
    else
        gA_dt      = A_dt
        gAt_dt     = At_dt
        gC_no_exp_ = C_no_exp_
        x_buf      = zeros(ComplexF64,n)
        vout_buf   = zeros(ComplexF64,n)
        term_buf   = zeros(ComplexF64,n)
        x_cpu_buf  = zeros(ComplexF64,n)
    end

    # Define the matrix-vector product function for eigs (no allocation inside)
    function matvec(x)
        return CONSTRUCT_C(x, gA_dt, gAt_dt, gC_no_exp_, n_t, freestream_disturbances, N_set_0, use_gpu, x_buf, vout_buf, term_buf, x_cpu_buf)
    end

    time_eigs = time()
    try
        (eigenvalues, eigenvectors) = eigs(LinearMap{ComplexF64}(matvec, n; ishermitian=true, isposdef=true), D_; nev=k, which=:LM, tol=1e-4, maxiter=10*k)
        elapsed_eigs = time() - time_eigs
        @printf("Eigenvalue solve time: %.6f seconds\n", elapsed_eigs)

        V = eigenvectors
        # The operator is theoretically symmetric so eigenvalues should be real,
        # but the Taylor-series matrix exponential introduces small imaginary
        # residuals.  Take the real part to avoid InexactError on ComplexF64→Float64.
        D = real.(eigenvalues)

        # Numerator gives energy of both real and complex parts of the mode. Each real and imaginary part have approximately same gain,
        # so divide by 2 to get the gain of each mode.
        D .= D ./ 2

        ## Print eigenvalues
        for i in 1:k
            @printf("eigenvalue %i = %f\n", i, real(D[i]))
        end
        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()  # Catch any async GPU errors
            CUDA.unsafe_free!(gA_dt)
            CUDA.unsafe_free!(gAt_dt)
            CUDA.unsafe_free!(gC_no_exp_)
            CUDA.unsafe_free!(x_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()  # Verify memory is deallocated
            CUDA.reclaim()
        end

        return (V, D, T_opt)
    catch e
        @warn "eigs did not converge: $e — returning NaN."
        V = fill(NaN, n, k)
        D = fill(NaN, k)

        elapsed = time() - time_eigs
        @printf("Elapsed time: %.6f seconds\n", elapsed)

        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()
            CUDA.unsafe_free!(gA_dt)
            CUDA.unsafe_free!(gAt_dt)
            CUDA.unsafe_free!(gC_no_exp_)
            CUDA.unsafe_free!(x_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()
            CUDA.reclaim()
        end

        return (V, D, T_opt)
    end
end
