# Tutorial: Non-Axisymmetric Stability of Pipe Flow (MATLAB)

Runtime: about 3 min (see the note on the device below)

## Goal

This tutorial shows how to analyse **non-axisymmetric disturbances** of an **axisymmetric base flow**. The base flow is laminar Hagen--Poiseuille flow in a pipe, and the disturbances have azimuthal wavenumber `m`:

$$
q'(x, r, \theta, t) = \hat q(x, r, t) \cos(m\theta), \qquad w'(x, r, \theta, t) = \hat w(x, r, t) \sin(m\theta),
$$

where $q' = (\rho', (\rho u)', (\rho v)', (\rho E)')$ and $w'$ is the azimuthal velocity. The base flow has no swirl, so this cos/sin pair decouples exactly and the linear operator stays **real**. The usual tools (eigenvalues, transient growth, receptivity) therefore run unchanged on it. The pipe is the simplest case with a symmetry axis and an analytic reference, which is why it is used to validate the implementation.

The workflow covers:
1. A **parametric sweep** over `m`, linearizing the same base flow once per wavenumber.
2. **Modal analysis**: least-stable eigenvalues for axial wavenumber $\alpha = 1$.
3. **Non-modal analysis**: optimal transient growth (lift-up) of streamwise-invariant disturbances, $\alpha = 0$, $m = 1$.
4. Comparison of every result with an independent spectral solution, `PIPE_FLOW_REFERENCE.m`.

---

## How azimuthal modes are switched on

Any `"3D-axisymmetric"` case accepts one more input:

```matlab
solution.stability_analysis.azimuthal_wavenumber = m;   % non-negative integer
[solution, L] = LINEARIZE_L(solution, chemistry);
```

When it is set, the linear state vector becomes `[rho; rho_u; rho_v; rho_E; rho_w; (shock)]`, with 5 fields per cell (`PERTURBATION_FIELDS`). For `m = 0` the swirl `rho_w` is decoupled, and the first four blocks reproduce the axisymmetric operator exactly. Leaving the field unset, or empty, gives the usual 4-variable operator.

The azimuthal terms are added by `FLUX_AZIMUTHAL_MODE` and `APPLY_BOUNDARY_CONDITIONS_AZIMUTHAL` (`utils/Operators/`):

- **Azimuthal fluxes:** their theta-derivatives in the mass, momentum and energy equations, including the m² viscous and heat-conduction terms, and the complete azimuthal-momentum equation in angular-momentum form.
- **Symmetry axis:** the ghost cell is the interior cell rotated by π, so every perturbation picks up a factor $(-1)^m$ relative to the axisymmetric mirror image.
- **Fitted shock:** its azimuthal tilt creates a post-shock swirl, $\hat w_s = \hat w_\infty - [u_n]\, (m / r_s)\, \hat\eta$.

Norms, transient growth and receptivity use the fifth variable automatically. The `w` kinetic energy is weighted like the radial velocity unless a fifth entry is given in `norms`.

---

## Step-by-Step Walkthrough

### Step 1 -- Inputs (`input_file.m`)

The geometry type `"pipe"` maps `chi` to the axial coordinate $x \in [0, L_x]$ and `eta` to the radius $r = y \in [0, R]$. The axis is at `eta0` (boundary condition `'symmetry'`) and the wall at `eta1`. The axial direction is periodic, so the domain holds axial wavenumbers $\alpha = 2\pi k / L_x$. With $L_x = 2\pi$ the least-stable modes of $\alpha = 1$ are resolved.

```matlab
solution.PDE_dimension = "3D-axisymmetric";
solution.curvilinear_mapping.boundary_type        = "pipe";
solution.curvilinear_mapping.Lx                   = 2*pi;
solution.curvilinear_mapping.R                    = 1;
solution.curvilinear_mapping.eta_refinement_power = 1.3;   % refinement towards the wall
solution.boundary_conditions.boundary_eta0.name = 'symmetry';           % axis
solution.boundary_conditions.boundary_eta1.name = 'no_slip_adiabatic';  % wall
solution.boundary_conditions.boundary_chi0.name = 'periodic';
solution.boundary_conditions.boundary_chi1.name = 'periodic';
solution.freestream.Mach = 0.05;   % centreline Mach number (nearly incompressible)
solution.freestream.Re   = 1000;   % Re = U_c R / nu
solution.stability_analysis.perturbation_magnitude = 1e-6;
solution.stability_analysis.azimuthal_wavenumber   = 1;
```

The base flow is the Hagen--Poiseuille profile $u = U_c (1 - r^2/R^2)$ with uniform density and pressure. As in the channel tutorial, it is linearized directly and no time marching is needed.

> **Perturbation magnitude.** At low Mach number the total energy scales like $1/M^2$. A finite-difference step of `1e-8` then loses digits in the pressure, which shows up in the slowest pressure-coupled modes (about 1% error at $M = 0.05$). `1e-6` is accurate to about `1e-4`.

### Step 2 -- Parametric sweep over `m` (`main.m`)

```matlab
for k = 1:numel(m_sweep)
    solution.stability_analysis.azimuthal_wavenumber = m_sweep(k);
    [solution_m, L] = LINEARIZE_L(solution, chemistry);
    reference = PIPE_FLOW_REFERENCE(Re, alpha, m_sweep(k), 60);
    modal(k).lambda_ref = conj(reference.lambda(1:n_modes));
    modal(k).lambda = arrayfun(@(lr) eigs(L, 1, lr), modal(k).lambda_ref);
end
```

Each wavenumber needs only a new `LINEARIZE_L`; the base flow is shared. The reference uses $\exp(i\alpha x)$. HYMOR's real operator contains both signs of $\alpha$, so the reference eigenvalues are searched as complex conjugates.

### Step 3 -- Transient growth, $\alpha = 0$, $m = 1$ (`main.m`)

```matlab
solution_tg.mesh.Nchi = 3;                  % only alpha = 0 resolved
solution_tg.mesh.Neta = 24;
solution_tg.curvilinear_mapping.Lx = 0.1;
solution_tg.freestream.Mach = 0.1;
solution_tg.stability_analysis.azimuthal_wavenumber = 1;
[solution_tg, L_tg] = LINEARIZE_L(solution_tg, chemistry);
[~, D_tg, T] = TRANSIENT_GROWTH_DOWNSTREAM(L_tg, solution_tg, 1, T);
```

- **Short period, three cells:** these resolve only streamwise-invariant disturbances. With $L_x = 2\pi$ and 3 cells, unresolved $\alpha \neq 0$ content raises the short-time gains by up to 10%.
- **Higher Mach number in this step:** the exponential time-stepper is acoustic, $\Delta t \propto M$. $M = 0.1$ halves its cost and changes $G$ by less than 0.2%.
- **Norm:** HYMOR maximizes Chu's energy. At these Mach numbers the optimal disturbance is purely kinetic, so $G$ compares directly with the incompressible kinetic-energy gain of the reference.

> **Device.** With `stability_analysis.device = "auto"` the exponential time-stepper runs on the GPU only for operators with at least $3\cdot10^5$ nonzeros: below that a GPU step is bound by kernel-launch latency. This operator has only a few hundred unknowns, so it runs on the CPU: 13 s and 66 s for $T = 10$ and $50$.

---

## Expected Results

Least-stable eigenvalues for $\alpha = 1$, $Re = 1000$, $M = 0.05$, on a 48 x 64 grid:

| m | HYMOR | Reference | Relative error |
|---|-------|-----------|----------------|
| 0 | -0.09025 + 0.90830i | -0.09035 + 0.91055i | 2.5e-3 |
| 0 | -0.15733 + 0.56455i | -0.15699 + 0.56615i | 2.8e-3 |
| 1 | -0.07086 + 0.84418i | -0.07086 + 0.84675i | 3.0e-3 |
| 1 | -0.09138 + 0.46759i | -0.09114 + 0.46914i | 3.3e-3 |
| 2 | -0.10426 + 0.80432i | -0.10434 + 0.80681i | 3.1e-3 |
| 2 | -0.13987 + 0.45175i | -0.13971 + 0.45312i | 2.9e-3 |
| 3 | -0.14412 + 0.76297i | -0.14429 + 0.76534i | 3.1e-3 |
| 3 | -0.16237 + 0.30317i | -0.16260 + 0.30400i | 2.5e-3 |

The m = 0 reference has two nearly degenerate modes near -0.0904 + 0.9106i, which this grid does not separate. The errors are second order: halving and doubling the grid gives about 1.2e-2 and 7e-4. At Re = 3000 the reference gives -0.041276 + 0.911466i for (m, alpha) = (1, 1), the classical value of Schmid & Henningson.

Optimal energy gain for $\alpha = 0$, $m = 1$ (3 x 24 grid, $M = 0.1$):

| T | HYMOR | Reference |
|---|-------|-----------|
| 10 | 16.90 | 16.96 |
| 50 | 72.01 | 72.34 |

The maximum over T, about 72 at $T \approx 50$, i.e. $G_{max}/Re^2 \approx 7.2 \times 10^{-5}$, is the lift-up optimum of pipe flow.

---

## The reference solution

`PIPE_FLOW_REFERENCE(Re, alpha, m, N)` solves the incompressible linearized Navier--Stokes equations about $U = 1 - r^2$ in primitive variables, using Chebyshev collocation on $r \in (0, 1]$. The parity of each variable across the axis (Trefethen, *Spectral Methods in MATLAB*) imposes regularity on the axis. It returns:

- the eigenvalues and velocity eigenvectors;
- `eval`, which interpolates an eigenvector at any radius;
- `gain`, the optimal kinetic-energy gain $G(T)$ from the eigenvector expansion (Reddy & Henningson).

It is independent of HYMOR's discretization. For $\alpha = 0$ its eigenvalues reduce to Bessel-function zeros:
- axial velocity: $-j_{m,n}^2 / Re$;
- cross-section Stokes modes: $-j_{m+1,n}^2 / Re$.

HYMOR's compressible operator also has the temperature (entropy) modes $-j'^2_{m,n} / (Re\, Pr)$ and the acoustic cross-section modes $\pm i\, j'_{m,n} / M$, which the incompressible reference does not contain. HYMOR reproduces both families to second order.
