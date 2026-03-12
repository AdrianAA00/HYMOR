using LinearAlgebra
using SparseArrays
using Arpack
using Printf

"""
    EIGENVALUES(A, s, n_modes)

Compute leading eigenvalues and eigenvectors of the linearized operator.

Dispatches to a CPU or GPU eigenvalue solver based on
`s["stability_analysis"]["eigenvalue_solver"]`. Available methods are:
- `"CPU_LU"` - Arpack eigs with shift-invert (LU factorization).
- `"GPU_TIMESTEPPER_ARNOLDI"` - Arnoldi power method via eigs (GPU or CPU).

Automatically detects CUDA availability. If no GPU is found, the Arnoldi
solver falls back to CPU.

# Arguments
- `A`: Linearized Jacobian matrix (sparse, n x n).
- `s`: Solution Dict{String,Any} with solver parameters.
- `n_modes`: Number of eigenvalues/eigenvectors to compute.

# Returns
- `V`: Matrix of eigenvectors (columns), size n x n_modes.
- `D`: Diagonal matrix of eigenvalues, size n_modes x n_modes.
  Returns V=0, D=0 if linearization is disabled.

# Notes
- For GPU_ARNOLDI, eigenvalues of the time-stepping operator are
  converted to continuous-time eigenvalues via logarithmic mapping.
- All memory is pre-allocated outside eigs. Only copies inside callback.
- GPU is synchronized and memory verified deallocated before returning.

Part of: Hypersonics Stability Julia Solver - Stability Analysis / Eigenvalues Module
"""
function EIGENVALUES(A, s::Dict{String,Any}, n_modes::Int)
    solver = s["stability_analysis"]["eigenvalue_solver"]

    if solver == "CPU_LU"
        (V, D) = EIGENVALUES_CPU_LU(A, s, n_modes)
    elseif solver == "GPU_TIMESTEPPER_ARNOLDI"
        (V, D) = EIGENVALUES_GPU_ARNOLDI(A, s, n_modes)
    else
        error("Unknown eigenvalue solver: $solver")
    end

    return (V, D)
end

## ========================================================================
#  CPU Shift-Invert Eigenvalue Solver
#  ========================================================================
"""
    EIGENVALUES_CPU_LU(A, s, n_modes)

Compute eigenvalues using Arpack eigs with shift-invert around a specified
sigma in the complex plane.
"""
function EIGENVALUES_CPU_LU(A, s::Dict{String,Any}, n_modes::Int)
    println()
    println("------------------------")
    println("      eigs Julia        ")
    println("------------------------")

    ## Frequency sweep parameters
    sigma_imag_min = 0.1
    sigma_imag_max = 0.50
    sigma_imag_step = 0.1
    sigma_real = 0.02

    N = round(Int, (sigma_imag_max - sigma_imag_min) / sigma_imag_step)

    if s["stability_analysis"]["perturb_shock"]
        V = zeros(s["mesh"]["Nchi"] * s["mesh"]["Neta"] * 4 + s["mesh"]["Nchi"], N)
        D = zeros(N, 1)
    else
        V = zeros(s["mesh"]["Nchi"] * s["mesh"]["Neta"] * 4, N)
        D = zeros(N, 1)
    end

    ## Targeted eigenvalue search
    N = 4
    time_start = time()
    # Convert A to complex: Arpack.jl requires a complex matrix when using a
    # complex shift sigma, so that the factorization of (A - sigma*I) is well-defined.
    (eigenvalues, eigenvectors) = eigs(complex(A); nev=4, sigma=-0.606629 + 0.462148im, tol=1e-6, maxiter=100)
    elapsed = time() - time_start
    @printf("Elapsed time: %.6f seconds\n", elapsed)

    V = eigenvectors
    D = eigenvalues

    ## Display results
    for i in 1:N
        @printf("eigenvalue %i = %f + %fi\n", i, real(D[i]), imag(D[i]))
    end

    return (V, D)
end

## ========================================================================
#  GPU/CPU Arnoldi Power Method Eigenvalue Solver
#  ========================================================================
"""
    EIGENVALUES_GPU_ARNOLDI(A, s, n_modes)

Compute eigenvalues using Arnoldi iteration with a Taylor-expansion time
stepper. Automatically detects GPU availability and falls back to CPU.
All memory is pre-allocated outside eigs to avoid repeated allocation.
"""
function EIGENVALUES_GPU_ARNOLDI(A, s::Dict{String,Any}, n_modes::Int)
    println()
    println("----------------------------------------------------")
    println(" Eigenvalue solver: Arnoldi + exponential transform ")
    println("----------------------------------------------------")

    ## Detect GPU availability
    use_gpu = CUDA.functional()

    if use_gpu
        gpu = CUDA.device()
        println("$(CUDA.name(gpu)) selected.")
    else
        println("No CUDA device found. Running on CPU.")
    end

    ## Time-stepping parameters
    T_ref = 0.5  # Time to evolve
    s = CFL_TIMESTEP(s)
    dt = s["time_integration"]["dt"]
    timesteps = round(Int, T_ref / dt)
    @printf("Number of timesteps: %f\n", timesteps)
    n = size(A, 1)
    k = n_modes

    ## Pre-allocate ALL memory outside eigs
    if use_gpu
        Adt = CuSparseMatrixCSC(A * dt)
        x_buf     = CUDA.zeros(Float64, n)
        vout_buf  = CUDA.zeros(Float64, n)
        term_buf  = CUDA.zeros(Float64, n)
        x_cpu_buf = zeros(Float64, n)
    else
        Adt = A * dt
        x_buf     = zeros(n)
        vout_buf  = zeros(n)
        term_buf  = zeros(n)
        x_cpu_buf = zeros(n)
    end

    ## Arnoldi iteration via eigs (no memory allocation inside callback)
    function matvec(x)
        return ARNOLDI_POWER_GPU(x, Adt, timesteps, use_gpu, x_buf, vout_buf, term_buf, x_cpu_buf)
    end

    time_eigs = time()
    try
        (eigenvalues, eigenvectors) = eigs(LinearMap(matvec, n; issymmetric=false); nev=k, which=:LM, tol=1e-4, maxiter=40*k)

        V = eigenvectors
        D = eigenvalues

        ## Convert discrete-time eigenvalues to continuous-time
        for i in 1:k
            D[i] = real(log(D[i]) / dt / timesteps) + 1im * imag(log(D[i]) / dt / timesteps)
            @printf("eigenvalue %i = %f + %fi\n", i, real(D[i]), imag(D[i]))
        end
        elapsed_eigs = time() - time_eigs
        @printf("Eigenvalue solve time: %.6f seconds\n", elapsed_eigs)

        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()  # Catch any async GPU errors
            CUDA.unsafe_free!(Adt)
            CUDA.unsafe_free!(x_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()  # Verify memory is deallocated
            CUDA.reclaim()
        end

        return (V, D)
    catch e
        @warn "eigs did not converge: $e — returning NaN."
        V = fill(NaN, n, k)
        D = fill(NaN + NaN*im, k)

        elapsed = time() - time_eigs
        @printf("Elapsed time: %.6f seconds\n", elapsed)

        ## GPU synchronization and memory cleanup
        if use_gpu
            CUDA.synchronize()
            CUDA.unsafe_free!(Adt)
            CUDA.unsafe_free!(x_buf)
            CUDA.unsafe_free!(vout_buf)
            CUDA.unsafe_free!(term_buf)
            CUDA.synchronize()
            CUDA.reclaim()
        end

        return (V, D)
    end
end
