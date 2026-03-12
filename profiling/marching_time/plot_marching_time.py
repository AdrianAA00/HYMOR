"""
plot_marching_time.py - Publication-quality performance comparison plot
                        for the time-marching solver (Julia vs MATLAB).

Reads timing results from timings_save/ and plots the cost per grid point
per time step per RK stage as a function of grid size.

Part of: Hypersonics Stability MATLAB Solver - Profiling Module
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_ITER = 1000       # number of time steps used in the benchmark
RK_STAGES = 4       # Explicit RK4
N_CORES = 8         # Number of CPU cores (threads)

DATA_DIR = Path(__file__).resolve().parent / "timings_save"
MATLAB_FILE = DATA_DIR / "timing_results_matlab.txt"
JULIA_FILE = DATA_DIR / "timing_results_julia.txt"
OUTPUT_FILE = Path(__file__).resolve().parent / "marching_time_performance.pdf"

# ---------------------------------------------------------------------------
# Helper: parse a timing results file
# ---------------------------------------------------------------------------
def parse_timing_file(filepath):
    """Return arrays of grid_sizes (Nchi*Neta), elapsed times (s), and labels."""
    grid_sizes = []
    elapsed = []
    labels = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("Grid"):
                continue
            parts = line.split()
            nchi, neta = parts[0].split("x")
            grid_sizes.append(int(nchi) * int(neta))
            elapsed.append(float(parts[-1]))
            labels.append(parts[0].replace("x", r"$\times$"))
    return np.array(grid_sizes), np.array(elapsed), labels

# ---------------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------------
grid_m, time_m, tick_labels = parse_timing_file(MATLAB_FILE)
grid_j, time_j, _ = parse_timing_file(JULIA_FILE)

# Normalise: seconds / grid_point / time_step / RK_stage
cost_m = time_m * N_CORES / (grid_m * N_ITER * RK_STAGES) * 1e6 # PID in microseconds
cost_j = time_j * N_CORES / (grid_j * N_ITER * RK_STAGES) * 1e6 # PID in microseconds

# ---------------------------------------------------------------------------
# Publication-quality plot
# ---------------------------------------------------------------------------
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

# Colors: MATLAB = blue, Julia = red  (consistent with eigen_solver_scaling)
color_matlab = "#1f77b4"
color_julia  = "#d62728"

# Markers: MATLAB = solid square, Julia = open circle
marker_matlab = "s"
marker_julia  = "o"

ls_matlab = "-"
ls_julia  = "--"

fig, ax = plt.subplots()

ax.plot(grid_m, cost_m, color=color_matlab, linestyle=ls_matlab, marker=marker_matlab,
        markerfacecolor=color_matlab, label="MATLAB")
ax.plot(grid_j, cost_j, color=color_julia, linestyle=ls_julia, marker=marker_julia,
        markerfacecolor="none", markeredgecolor=color_julia, label="Julia")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel(r"Grid size $N_\chi \times N_\eta$")
ax.set_ylabel(r"PID ($\mu$s / DOF)")

# Custom x-tick labels showing NxN dimensions
ax.set_xticks(grid_m)
ax.set_xticklabels(tick_labels, rotation=35, ha="right")
ax.xaxis.set_minor_locator(plt.NullLocator())

# Show readable numeric labels on y-axis (data spans < 1 decade)
from matplotlib.ticker import ScalarFormatter
yfmt = ScalarFormatter()
yfmt.set_scientific(True)
yfmt.set_powerlimits((-1, 1))
ax.yaxis.set_major_formatter(yfmt)
ax.yaxis.set_minor_formatter(yfmt)

ax.legend(frameon=True, fancybox=False, edgecolor="black")
ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.7)
ax.tick_params(direction="in", which="both", top=True, right=True)

fig.tight_layout()
fig.savefig(OUTPUT_FILE)
fig.savefig(OUTPUT_FILE.with_suffix(".png"))
print(f"Saved  {OUTPUT_FILE}")
print(f"Saved  {OUTPUT_FILE.with_suffix('.png')}")
plt.show()
