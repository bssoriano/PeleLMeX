#!/usr/bin/env python3
"""
Mid-plane (z=0) scatter plot of Temperature vs. Progress Variable, colored by
the normalized log scalar dissipation rate, with the fully burning (low-strain)
flamelet solution overlaid as a red dashed ideal-combustion reference --
reproducing the Fig. 12 style from the reference paper.

Single snapshot:
    python3 plot_scatter_T_C_midplane.py [pltdir]

Time-ensemble (pools points from many snapshots into one, fuller scatter cloud):
    python3 plot_scatter_T_C_midplane.py --ensemble
    python3 plot_scatter_T_C_midplane.py --ensemble --start 112707 --end 133981 --stride 2
"""
import os
import re
import warnings
import argparse

import yt
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Data extraction ---
T_FIELD = ("boxlib", "temp")
C_FIELD = ("boxlib", "Y(PROG)")
CHI_FIELD = ("boxlib", "chi_sgs")   # solver's stored field is RHO*chi_sgs, not chi_sgs -- see extract_midplane_data
DENSITY_FIELD = ("boxlib", "density")
VOLFRAC_FIELD = ("boxlib", "volFrac")
BUFF_SIZE = (400, 1200)     # mid-plane resampling resolution (dx == dy)
VOLFRAC_MIN = 0.5           # exclude EB (nozzle) solid/cut cells
ENSEMBLE_MIN_LEVEL = 3      # require this many AMR levels so pooled snapshots share resolution

# --- Flamelet reference (streams match input.3d-regt; "O2:1.0, N2:3.76" is a
# --- mole ratio, so streams are set via TPX, not TPY) ---
MECHANISM = os.path.join(SCRIPT_DIR, "generate_manifold", "drm19.yaml")
PRESSURE = 101325.0
FUEL_T, FUEL_COMP = 298.0, "CH4:1"
OX_T, OX_COMP = 298.0, "O2:1.0, N2:3.76"
FLAMELET_WIDTH = 0.6
FLAMELET_MDOT = 0.1
FLAMELET_TRANSPORT = "unity-Lewis-number"
PROGRESS_SPECIES = ["CO2", "H2O", "CO", "H2"]   # C = Y_CO2 + Y_H2O + Y_CO + Y_H2

# --- Rendering ---
MAX_SCATTER_POINTS = 400000
SEED = 0
SCATTER_POINT_SIZE = 12

XLIM = (0.0, 0.3)
YLIM = (298.0, 2200.0)


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
    C = np.array(frb[C_FIELD]).ravel()
    rho_chi = np.array(frb[CHI_FIELD]).ravel()
    rho = np.array(frb[DENSITY_FIELD]).ravel()
    volfrac = np.array(frb[VOLFRAC_FIELD]).ravel()

    fluid = volfrac >= volfrac_min
    C = np.clip(C[fluid], 0.0, None)
    T = T[fluid]
    # PeleLMeX's "chi_sgs" plotfile variable is actually rho*chi_sgs (see
    # PeleLMeX_Forces.cpp's linear-relaxation model: "rho chi_sgs = C_chi *
    # mu_t / Delta^2 * Variance", needed as-is to source the conserved
    # rho*Z''^2 variance equation) -- divide by density to get chi_sgs itself.
    chi = np.clip(rho_chi[fluid] / rho[fluid], 0.0, None)
    return C, T, chi


def collect_ensemble_data(plotfiles):
    C_all, T_all, chi_all = [], [], []
    for i, pf in enumerate(plotfiles):
        name = os.path.basename(pf)
        try:
            C, T, chi = extract_midplane_data(pf)
        except Exception as e:
            print(f"  [{i + 1}/{len(plotfiles)}] {name}: skipped ({e})")
            continue
        print(f"  [{i + 1}/{len(plotfiles)}] {name}: {C.size} fluid points")
        C_all.append(C)
        T_all.append(T)
        chi_all.append(chi)
    return np.concatenate(C_all), np.concatenate(T_all), np.concatenate(chi_all)


def compute_flamelet_reference(mechanism, fuel_T, fuel_comp, ox_T, ox_comp, pressure,
                                width=FLAMELET_WIDTH, mdot=FLAMELET_MDOT,
                                transport=FLAMELET_TRANSPORT):
    """Solve a single low-strain counterflow diffusion flame and return its
    (C, T) trajectory ordered by mixture fraction Z, *not* sorted by C.

    C(Z) is unimodal (rises from the oxidizer side, peaks near stoichiometric,
    falls back toward zero on the fuel-rich side), so a given C generally
    corresponds to two different T values -- one per side of the peak.
    Preserving the Z-ordering lets the plotted line fold back on itself
    correctly instead of a sort-by-C artifact connecting unrelated branches.
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
    species_idx = [flame.gas.species_names.index(s) for s in PROGRESS_SPECIES]
    C = flame.Y[species_idx, :].sum(axis=0)

    order = np.argsort(Z)
    return C[order], T[order]


def normalize_log(values, calibration_mask=None, low_percentile=1.0, high_percentile=99.0):
    """log10-scale and normalize to [0, 1] (1 = highest, 0 = lowest), matching
    the convention used for point density in Fig. 10/14.

    Most of the mid-plane is unreacted (C=0, chi_sgs~0 or below float
    precision), so a plain min/max would let a handful of near-zero cells
    set the floor and compress the whole active-flame range into a sliver
    near 1. Calibrating the [low_percentile, high_percentile] range on just
    the reacting subset (calibration_mask) keeps the color scale informative
    there; values are clipped to [0, 1] afterward so the unreacted background
    just floors out at 0 instead of skewing the scale.
    """
    positive = values[values > 0]
    floor = positive.min() if positive.size else 1e-30
    log_values = np.log10(np.maximum(values, floor))
    calibration = log_values if calibration_mask is None else log_values[calibration_mask]
    lo, hi = np.percentile(calibration, [low_percentile, high_percentile])
    return np.clip((log_values - lo) / (hi - lo), 0.0, 1.0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Mid-plane T-C scatter plot colored by normalized log scalar "
                    "dissipation rate, with a fully burning flamelet reference curve.")
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
    parser.add_argument("-o", "--output", default=None, help="Output PDF path")
    args = parser.parse_args()

    if args.ensemble:
        plotfiles = find_plotfiles(SCRIPT_DIR, min_level=args.min_level,
                                    start=args.start, end=args.end, stride=args.stride)
        if not plotfiles:
            raise FileNotFoundError("No plotfiles matched the --ensemble selection.")
        print(f"Pooling {len(plotfiles)} snapshots "
              f"({os.path.basename(plotfiles[0])} .. {os.path.basename(plotfiles[-1])}):")
        C, T, chi = collect_ensemble_data(plotfiles)
        tag = f"ensemble_{os.path.basename(plotfiles[0])}-{os.path.basename(plotfiles[-1])}_n{len(plotfiles)}"
    else:
        plotfile = args.plotfile or find_plotfiles(SCRIPT_DIR)[-1]
        print(f"Loading {plotfile} ...")
        C, T, chi = extract_midplane_data(plotfile)
        tag = os.path.basename(os.path.normpath(plotfile))

    print(f"Total {C.size} fluid-phase points from the mid-plane (z=0).")

    chi_norm = normalize_log(chi, calibration_mask=(C > 1e-4))

    rng = np.random.default_rng(SEED)
    if C.size > MAX_SCATTER_POINTS:
        idx = rng.choice(C.size, size=MAX_SCATTER_POINTS, replace=False)
        C, T, chi_norm = C[idx], T[idx], chi_norm[idx]

    order = np.argsort(chi_norm)
    C, T, chi_norm = C[order], T[order], chi_norm[order]

    print("Solving the fully burning (low-strain) flamelet reference...")
    C_fl, T_fl = compute_flamelet_reference(MECHANISM, FUEL_T, FUEL_COMP,
                                             OX_T, OX_COMP, PRESSURE)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sc = ax.scatter(C, T, c=chi_norm, cmap="viridis", vmin=0.0, vmax=1.0, s=SCATTER_POINT_SIZE,
                    alpha=0.7, linewidths=0, rasterized=True)
    ax.plot(C_fl, T_fl, "r--", linewidth=1.8,
            label="fully burning flamelet\n(cantera)")

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_xlabel(r"Progress Variable, $C$")
    ax.set_ylabel(r"Temperature [K]")
    ax.legend(loc="upper right", fontsize=8, frameon=False)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02, ticks=[0, 0.25, 0.5, 0.75, 1])
    cbar.set_label(r"Normalized log scalar dissipation rate, $\chi_{sgs}$")

    fig.tight_layout()
    fig_dir = os.path.join(SCRIPT_DIR, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    out_path = args.output or os.path.join(fig_dir, f"scatter_T_C_midplane_{tag}.pdf")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")
    plt.show()
