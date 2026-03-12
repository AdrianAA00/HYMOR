import numpy as np
import matplotlib.pyplot as plt
import os

# --- Data loading ---
timings_dir = os.path.join(os.path.dirname(__file__), "timings_save")

def load_timing(filename):
    """Load timing file and return grid_size and time arrays."""
    filepath = os.path.join(timings_dir, filename)
    grid_size = []
    time = []
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.split()
            grid_size.append(int(parts[2]))
            time.append(float(parts[4]))
    return np.array(grid_size), np.array(time)

# Load all data
eig_grid_m,   eig_time_m   = load_timing("eigenvalues_timings_MATLAB.txt")
eig_grid_j,   eig_time_j   = load_timing("eigenvalues_timings_julia.txt")
tg_grid_m,    tg_time_m    = load_timing("transient_growth_downstream_timings_MATLAB.txt")
tg_grid_j,    tg_time_j    = load_timing("transient_growth_downstream_timings_julia.txt")
fs_grid_m,    fs_time_m    = load_timing("freestream_receptivity_timings_MATLAB.txt")
fs_grid_j,    fs_time_j    = load_timing("freestream_receptivity_timings_julia.txt")

# --- Publication formatting ---
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 11,
    "axes.labelsize": 12,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
    "figure.figsize": (5.5, 4.5),
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# Colors: same solver type = same color
color_eig = "#1f77b4"  # blue
color_tg  = "#d62728"  # red
color_fs  = "#2ca02c"  # green

# Markers: MATLAB = solid markers, Julia = open markers
marker_matlab = "s"  # square
marker_julia  = "o"  # circle

ls_matlab = "-"
ls_julia  = "--"

fig, ax = plt.subplots()

# Eigenvalues
ax.plot(eig_grid_m, eig_time_m, color=color_eig, linestyle=ls_matlab, marker=marker_matlab,
        markerfacecolor=color_eig, label="Modal analysis (MATLAB)")
ax.plot(eig_grid_j, eig_time_j, color=color_eig, linestyle=ls_julia, marker=marker_julia,
        markerfacecolor="none", markeredgecolor=color_eig, label="Modal analysis (Julia)")

# Transient growth
ax.plot(tg_grid_m, tg_time_m, color=color_tg, linestyle=ls_matlab, marker=marker_matlab,
        markerfacecolor=color_tg, label="Transient growth (MATLAB)")
ax.plot(tg_grid_j, tg_time_j, color=color_tg, linestyle=ls_julia, marker=marker_julia,
        markerfacecolor="none", markeredgecolor=color_tg, label="Transient growth (Julia)")

# Freestream receptivity
ax.plot(fs_grid_m, fs_time_m, color=color_fs, linestyle=ls_matlab, marker=marker_matlab,
        markerfacecolor=color_fs, label="Freestream receptivity (MATLAB)")
ax.plot(fs_grid_j, fs_time_j, color=color_fs, linestyle=ls_julia, marker=marker_julia,
        markerfacecolor="none", markeredgecolor=color_fs, label="Freestream receptivity (Julia)")

# --- Reference power-law line: t ∝ N^1.5 ---
# Anchor: fit constant C so the line passes through the geometric mean
# of all datasets at their median grid point, keeping it visually central.
exponent = 1.5
all_grids = np.concatenate([eig_grid_m, eig_grid_j,
                             tg_grid_m,  tg_grid_j,
                             fs_grid_m,  fs_grid_j])
all_times = np.concatenate([eig_time_m, eig_time_j,
                             tg_time_m,  tg_time_j,
                             fs_time_m,  fs_time_j])

# Least-squares anchor in log-space: minimise sum of (log t - log C - 1.5 log N)^2
# => log C = mean(log t - 1.5 * log N), then scale down one decade for readability
log_C = np.mean(np.log(all_times) - exponent * np.log(all_grids))
C_anchor = np.exp(log_C) * 3   # shift below the data cloud so it doesn't occlude

ref_N = np.logspace(np.log10(all_grids.min()),
                    np.log10(all_grids.max()), 200)
ref_t = C_anchor * ref_N ** exponent

ax.plot(ref_N, ref_t,
        color="0.35", linestyle=(0, (5, 4, 1, 4)),   # dash-dot: 5 on, 4 off, 1 on, 4 off
        linewidth=1, zorder=1)

# Slope annotation: place at the right-hand end of the reference line
x_ann  = ref_N[-1] * 0.03          # a bit left of the right edge
y_ann  = (C_anchor * 1.2) * x_ann ** exponent
# Angle of the line in display coordinates — computed after draw so we use a
# fixed rotation consistent with log-log slope-triangle convention instead.
ax.annotate(r"$\propto (N_\chi \times N_\eta)^{1.5}$",
            xy=(x_ann, y_ann),
            xytext=(-14, 6),   # negative x offset nudges label left
            textcoords="offset points",
            fontsize=12, color="0.35",
            ha="center", va="bottom",
            rotation=0)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"Grid size $N_\chi \times N_\eta$")
ax.set_ylabel(r"Time [s]")
ax.legend(frameon=True, fancybox=False, edgecolor="black")
ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.7)
ax.tick_params(direction="in", which="both", top=True, right=True)

fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(__file__), "eigen_solver_scaling.pdf"))
fig.savefig(os.path.join(os.path.dirname(__file__), "eigen_solver_scaling.png"))
plt.show()