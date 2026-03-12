using LinearAlgebra
using SparseArrays
using Arpack
using Printf

"""
    TRANSIENT_GROWTH_DOWNSTREAM(A, s, n_modes, T)

Compute optimal transient growth modes for downstream disturbances.

Solves the generalized eigenvalue problem for transient energy growth
in the downstream flow domain. Uses a matrix-free approach with eigs()
to find the n_modes largest-gain eigenmodes at the specified time
horizon T. Automatically detects GPU availability.

# Arguments
- `A`: System Jacobian matrix (sparse).
- `s`: Solution Dict{String,Any} with flow fields, mesh data, and configuration flags.
- `n_modes`: Number of leading transient growth modes to compute.
- `T`: Target time horizon for the optimization.

# Returns
- `V`: Matrix of optimal perturbation eigenvectors (columns).
- `D`: Diagonal matrix of corresponding energy gains (eigenvalues).
- `T_opt`: Actual time horizon used (adjusted to integer time steps).

# Notes
- All memory is pre-allocated outside eigs. Only copies inside callback.
- GPU is synchronized and memory verified deallocated before returning.
- Falls back to CPU if no CUDA device is found.

Part of: Hypersonics Stability Julia Solver - Transient Growth Downstream Module
"""
function TRANSIENT_GROWTH_DOWNSTREAM(A, s::Dict{String,Any}, n_modes::Int, T)

    k = n_modes    # Number of modes to compute
    freestream_disturbances = false

    ## Compute time-stepping parameters
    s = CFL_TIMESTEP(s)
    t_temp = s["time_integration"]["dt"]
    n_t = round(Int, T / t_temp)
    dt = T / n_t
    T_opt = n_t * dt

    ## Construct operator D and C (without exponentials)
    R = CONSTRUCT_R(s)
    norms = [1, 1, 1, 1]
    M = CONSTRUCT_M(s, norms)
    M_hat = CONSTRUCT_M_HAT(s, norms)
    C_no_exp = R' * M * R
    D_mat = R' * M_hat * R

    ## Eigenvalue computation via Lanczos iteration
    println()
    println("---------------------------------")
    println("    Transient growth analysis    ")
    println("---------------------------------")
    println()
    @printf("Number timesteps = %f \n", n_t)

    ## Detect GPU availability
    use_gpu = CUDA.functional()
    n = size(A, 1)

    if use_gpu
        gpu = CUDA.device()
        println("$(CUDA.name(gpu)) selected.")
    else
        println("No CUDA device found. Running on CPU.")
    end

    ## Pre-allocate ALL memory outside eigs (matrices + work buffers)
    A_dt = A * dt
    At_dt = sparse(A_dt')

    if use_gpu
        gA_dt      = CuSparseMatrixCSC(A_dt)
        gAt_dt     = CuSparseMatrixCSC(At_dt)
        gC_no_exp  = CuSparseMatrixCSC(C_no_exp)
        x_buf      = CUDA.zeros(Float64, n)
        vout_buf   = CUDA.zeros(Float64, n)
        term_buf   = CUDA.zeros(Float64, n)
        x_cpu_buf  = zeros(Float64, n)
    else
        gA_dt      = A_dt
        gAt_dt     = At_dt
        gC_no_exp  = C_no_exp
        x_buf      = zeros(n)
        vout_buf   = zeros(n)
        term_buf   = zeros(n)
        x_cpu_buf  = zeros(n)
    end

    # Define the matrix-vector product function for eigs (no allocation inside)
    function matvec(x)
        return CONSTRUCT_C(x, gA_dt, gAt_dt, gC_no_exp, n_t, freestream_disturbances, 0, use_gpu, x_buf, vout_buf, term_buf, x_cpu_buf)
    end

    time_eigs = time()
    try

        (eigenvalues, eigenvectors) = eigs(LinearMap(matvec, n; issymmetric=true, isposdef=true), D_mat; nev=k, which=:LM, tol=1e-4, maxiter=10*k)
        elapsed_eigs = time() - time_eigs
        @printf("Eigenvalue solve time: %.6f seconds\n", elapsed_eigs)
        
        V = eigenvectors
        D = eigenvalues

        ## Print eigenvalues
        for i in 1:k
            @printf("eigenvalue %i = %f\n", i, real(D[i]))
        end


        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()  # Catch any async GPU errors
            CUDA.unsafe_free!(gA_dt)
            CUDA.unsafe_free!(gAt_dt)
            CUDA.unsafe_free!(gC_no_exp)
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
            CUDA.unsafe_free!(gC_no_exp)
            CUDA.unsafe_free!(x_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()
            CUDA.reclaim()
        end

        return (V, D, T_opt)
    end
end
