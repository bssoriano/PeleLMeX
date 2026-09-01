#!/usr/bin/env python3
"""
Mid-plane (z=0) scatter plot of Temperature vs. Mixture Fraction, colored by
normalized log point density, with the fully burning (low-strain) flamelet
solution overlaid as a red dashed ideal-combustion reference -- reproducing
the Fig. 10 style from the reference paper.

Single snapshot:
    python3 plot_scatter_T_Z_midplane.py [pltdir]

Time-ensemble (pools points from many snapshots into one, fuller scatter cloud):
    python3 plot_scatter_T_Z_midplane.py --ensemble
    python3 plot_scatter_T_Z_midplane.py --ensemble --start 112707 --end 133981 --stride 2
"""
import os
import re
import warnings
import argparse

import yt
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Data extraction ---
T_FIELD = ("boxlib", "temp")
Z_FIELD = ("boxlib", "Y(ZMIX)")
VOLFRAC_FIELD = ("boxlib", "volFrac")
BUFF_SIZE = (400, 1200)     # mid-plane resampling resolution (dx == dy)
VOLFRAC_MIN = 0.5           # exclude EB (nozzle) solid/cut cells
ENSEMBLE_MIN_LEVEL = 3      # require this many AMR levels so pooled snapshots share resolution

# --- Flamelet reference (streams match input.3d-regt; "O2:1.0, N2:3.76" is a
# --- mole ratio, so streams are set via TPX, not TPY) ---
MECHANISM = os.path.join(SCRIPT_DIR, "generate_manifold", "drm19.yaml")
PRESSURE = 101325.0         # prob.P_mean [Pa]
FUEL_T, FUEL_COMP = 298.0, "CH4:1"                 # prob.jet_T, prob.jet_Yfuel=1.0
OX_T, OX_COMP = 298.0, "O2:1.0, N2:3.76"           # prob.T_mean (air)
FLAMELET_WIDTH = 0.6        # matches generate_manifold/get_nonpremixed_flamelets.py
FLAMELET_MDOT = 0.1         # kg/m2/s on both inlets -- lowest-strain (most fully-burning) case
FLAMELET_TRANSPORT = "unity-Lewis-number"

# --- Density estimate & rendering ---
DENSITY_BINS = 200
DENSITY_SMOOTH_SIGMA = 1.5
MAX_SCATTER_POINTS = 400000  # most of the mid-plane is unmixed air; keep the rare high-Z tail intact
SEED = 0
SCATTER_POINT_SIZE = 22     # marker area in points^2 (matplotlib's `s`)

XLIM = (0.0, 1.0)
YLIM = (0.0, 2000.0)


def find_plotfiles(base_dir, min_level=None, start=None, end=None, stride=1):
    """Sorted plt##### paths in base_dir, optionally filtered by AMR depth
    (min_level) and bounded/strided by their numeric step index."""
    candidates = []
    for name in os.listdir(base_dir):
        m = re.fullmatch(r"plt(\d+)", name)
        if not m or not os.path.isdir(os.path.join(base_dir, name)):
            continue
        num = int(m.group(1))
        if start is not None and num < start:
            continue
        if end is not None and num > end:
            continue
        if min_level is not None and not os.path.isdir(os.path.join(base_dir, name, f"Level_{min_level}")):
            continue
        candidates.append((num, name))
    candidates.sort()
    return [os.path.join(base_dir, name) for _, name in candidates[::stride]]


def extract_midplane_data(plotfile, buff_size=BUFF_SIZE, volfrac_min=VOLFRAC_MIN):
    ds = yt.load(plotfile)
    center = [((ds.domain_left_edge[0] + ds.domain_right_edge[0]) / 2).v,
              ((ds.domain_left_edge[1] + ds.domain_right_edge[1]) / 2).v,
              0.0]

    slc = yt.SlicePlot(ds, "z", T_FIELD, center=center)
    slc.set_buff_size(buff_size)
    frb = slc.frb

    T = np.array(frb[T_FIELD]).ravel()
    Z = np.array(frb[Z_FIELD]).ravel()
    volfrac = np.array(frb[VOLFRAC_FIELD]).ravel()

    fluid = volfrac >= volfrac_min
    Z = np.clip(Z[fluid], 0.0, 1.0)
    T = T[fluid]
    return Z, T


def collect_ensemble_data(plotfiles):
    Z_all, T_all = [], []
    for i, pf in enumerate(plotfiles):
        name = os.path.basename(pf)
        try:
            Z, T = extract_midplane_data(pf)
        except Exception as e:
            print(f"  [{i + 1}/{len(plotfiles)}] {name}: skipped ({e})")
            continue
        print(f"  [{i + 1}/{len(plotfiles)}] {name}: {Z.size} fluid points")
        Z_all.append(Z)
        T_all.append(T)
    return np.concatenate(Z_all), np.concatenate(T_all)


def compute_flamelet_reference(mechanism, fuel_T, fuel_comp, ox_T, ox_comp, pressure,
                                width=FLAMELET_WIDTH, mdot=FLAMELET_MDOT,
                                transport=FLAMELET_TRANSPORT):
    """Solve a single low-strain counterflow diffusion flame: the "fully
    burning flamelet solution" used as the ideal-combustion reference.

    This is deliberately a real flame solve rather than a per-Z chemical
    equilibrium sweep at the mixing enthalpy h_mix(Z): checked numerically,
    this flame's actual local enthalpy departs from linear mixing (a lean-side
    deficit / rich-side surplus of up to ~30000 J/kg for this fuel/oxidizer
    pair), so a naive per-Z equilibrium at h_mix(Z) is *not* an upper bound on
    it -- the flamelet exceeded such a naive reference by up to 578 K on the
    rich side in testing. A real flame solve sidesteps that entirely.
    """
    fuel_gas = ct.Solution(mechanism)
    fuel_gas.TPX = fuel_T, pressure, fuel_comp
    ox_gas = ct.Solution(mechanism)
    ox_gas.TPX = ox_T, pressure, ox_comp

    gas = ct.Solution(mechanism)
    flame = ct.CounterflowDiffusionFlame(gas, width=width)
    flame.P = pressure
    flame.fuel_inlet.Y = fuel_gas.Y
    flame.fuel_inlet.T = fuel_T
    flame.fuel_inlet.mdot = mdot
    flame.oxidizer_inlet.Y = ox_gas.Y
    flame.oxidizer_inlet.T = ox_T
    flame.oxidizer_inlet.mdot = mdot
    flame.set_refine_criteria(ratio=3.0, slope=0.15, curve=0.15, prune=0.05)
    flame.flame.set_steady_tolerances(default=[1e-6, 1e-12])
    flame.max_grid_points = 8192
    flame.transport_model = transport

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        flame.solve(loglevel=0, auto=True)

    Z = flame.mixture_fraction("N")
    T = flame.T
    order = np.argsort(Z)
    return Z[order], T[order]


def compute_point_density(x, y, bins=DENSITY_BINS, smooth_sigma=DENSITY_SMOOTH_SIGMA):
    """Local population density per point, log10-scaled and min-max
    normalized to [0, 1] (1 = densest clustering, 0 = sparsest), matching the
    convention in the reference paper's Figure 10."""
    H, xedges, yedges = np.histogram2d(x, y, bins=bins)
    H = gaussian_filter(H, sigma=smooth_sigma)

    ix = np.clip(np.digitize(x, xedges) - 1, 0, H.shape[0] - 1)
    iy = np.clip(np.digitize(y, yedges) - 1, 0, H.shape[1] - 1)
    density = H[ix, iy]

    log_density = np.log10(np.maximum(density, density[density > 0].min()))
    return (log_density - log_density.min()) / (log_density.max() - log_density.min())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Mid-plane T-Z scatter plot colored by normalized log point "
                    "density, with a fully burning flamelet reference curve.")
    parser.add_argument("plotfile", nargs="?", default=None,
                        help="Single PeleLMeX plt directory (default: latest plt##### found). "
                             "Ignored if --ensemble is set.")
    parser.add_argument("--ensemble", action="store_true",
                        help="Pool mid-plane points from many snapshots instead of just one.")
    parser.add_argument("--min-level", type=int, default=ENSEMBLE_MIN_LEVEL,
                        help="Ensemble mode: only use plotfiles with at least this many "
                             "AMR levels, so pooled snapshots share resolution (default: %(default)s).")
    parser.add_argument("--start", type=int, default=None, help="Ensemble mode: minimum plt step number.")
    parser.add_argument("--end", type=int, default=None, help="Ensemble mode: maximum plt step number.")
    parser.add_argument("--stride", type=int, default=1,
                        help="Ensemble mode: use every Nth matching snapshot (default: %(default)s).")
    parser.add_argument("-o", "--output", default=None, help="Output PNG path")
    args = parser.parse_args()

    if args.ensemble:
        plotfiles = find_plotfiles(SCRIPT_DIR, min_level=args.min_level,
                                    start=args.start, end=args.end, stride=args.stride)
        if not plotfiles:
            raise FileNotFoundError("No plotfiles matched the --ensemble selection.")
        print(f"Pooling {len(plotfiles)} snapshots "
              f"({os.path.basename(plotfiles[0])} .. {os.path.basename(plotfiles[-1])}):")
        Z, T = collect_ensemble_data(plotfiles)
        tag = f"ensemble_{os.path.basename(plotfiles[0])}-{os.path.basename(plotfiles[-1])}_n{len(plotfiles)}"
    else:
        plotfile = args.plotfile or find_plotfiles(SCRIPT_DIR)[-1]
        print(f"Loading {plotfile} ...")
        Z, T = extract_midplane_data(plotfile)
        tag = os.path.basename(os.path.normpath(plotfile))

    print(f"Total {Z.size} fluid-phase points from the mid-plane (z=0).")

    density = compute_point_density(Z, T)

    rng = np.random.default_rng(SEED)
    if Z.size > MAX_SCATTER_POINTS:
        idx = rng.choice(Z.size, size=MAX_SCATTER_POINTS, replace=False)
        Z, T, density = Z[idx], T[idx], density[idx]

    order = np.argsort(density)
    Z, T, density = Z[order], T[order], density[order]

    print("Solving the fully burning (low-strain) flamelet reference...")
    Z_fl, T_fl = compute_flamelet_reference(MECHANISM, FUEL_T, FUEL_COMP,
                                             OX_T, OX_COMP, PRESSURE)

    fig, ax = plt.subplots(figsize=(8, 7.2))
    sc = ax.scatter(Z, T, c=density, cmap="BuPu", vmin=0.0, vmax=1.0, s=SCATTER_POINT_SIZE, alpha=0.7,
                    linewidths=0, rasterized=True)
    ax.plot(Z_fl, T_fl, "r--", linewidth=1.8,
            label="fully burning flamelet\n(cantera)")

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([298, 1000, 2300])
    ax.set_xlabel(r"Mixture Fraction, $Z$")
    ax.set_ylabel(r"Temperature [K]")
    ax.legend(loc="upper right", fontsize=8, frameon=False)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02, ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label("Normalized log point density")

    fig.tight_layout()
    fig_dir = os.path.join(SCRIPT_DIR, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    out_path = args.output or os.path.join(fig_dir, f"scatter_T_Z_midplane_{tag}.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()
