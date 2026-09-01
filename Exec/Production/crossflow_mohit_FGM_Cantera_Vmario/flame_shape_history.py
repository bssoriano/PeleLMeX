#!/usr/bin/env python3
"""
Flame shape (length & width) history for the crossflow FGM case.

Method: threshold the temperature field at T = 800 K, following Mohit
et al. (2026, Flow Turbul. Combust. 116:6), Sec. 3.1, after Bradley et al.
(2016):
  - "length" = farthest downstream extent of the T=800 K envelope,
               measured from the nozzle location.
  - "width"  = diameter of the cylinder (axis along the downstream
               direction) that encloses the entire T=800 K envelope, i.e.
               the diameter of the envelope's footprint projected onto
               the plane perpendicular to the downstream direction.
    (Implemented as the point-set diameter -- max pairwise distance over
    the convex hull -- of that footprint; a robust, dependency-light
    stand-in for the minimal-enclosing-circle diameter the paper uses.)

The paper's own headline numbers (its Table 2/3) are this measurement
applied to the TIME-AVERAGED temperature field over a statistically
steady window. This script reproduces that ("TIME-AVERAGED FLAME SHAPE"
below) and, in addition, applies the same measurement to every individual
snapshot so you can see the flame-shape time history/fluctuations.

Axis convention for THIS case (see input.3d-regt / pelelmex_prob.H):
  y = downstream / crossflow direction (Inflow at y=0, Outflow at y=ymax)
  x = jet direction (fuel exits the EB nozzle blowing along +x)
  z = spanwise direction
All lengths in meters, times in seconds.

Usage:
    python3 flame_shape_history.py [first_plt] [last_plt]
Both arguments are optional; if omitted, FIRST_PLT/LAST_PLT below are
used. Every plt### directory in the current folder whose step number
falls in [first, last] (inclusive) is processed.
"""
import matplotlib
matplotlib.use("Agg")

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist
from scipy.ndimage import label
import yt

yt.funcs.mylog.setLevel(50)  # silence yt INFO spam

# ---------------------------------------------------------------------
# CONFIGURATION - edit these, or override first/last plt on the command line
# ---------------------------------------------------------------------
FIRST_PLT = "plt18520"     # first plotfile of the averaging/analysis window
LAST_PLT =  "plt1011701"      # last plotfile of the averaging/analysis window (edit as the run progresses)

T_FLAME = 800.0            # K, flame-boundary threshold (Bradley et al. 2016)
JET_Y = 0.762               # m, downstream location of the nozzle (10 D) -> length origin
D_NOZZLE = 2.0 * 0.0381     # m, nozzle diameter (2 * prob.jet_radius, input.3d-regt)

GRID_LEVEL = 2              # AMR level every snapshot is resampled onto (0=coarsest)

OUT_DIR = "flame_shape_results"
# ---------------------------------------------------------------------


def get_plotfiles_in_range(first_name, last_name):
    start_num = int(re.search(r"\d+", first_name).group())
    end_num = int(re.search(r"\d+", last_name).group())
    all_dirs = [d for d in os.listdir(".") if os.path.isdir(d) and d.startswith("plt")]
    selected = []
    for d in all_dirs:
        m = re.search(r"\d+", d)
        if m and start_num <= int(m.group()) <= end_num:
            selected.append(d)
    return sorted(selected, key=lambda d: int(re.search(r"\d+", d).group()))


def enclosing_diameter(px, pz, max_points=4000):
    """Point-set diameter (max pairwise distance) of a 2D point cloud,
    computed via its convex hull for speed. Used as the flame 'width'."""
    pts = np.unique(np.column_stack([px, pz]), axis=0)
    if len(pts) < 2:
        return 0.0
    try:
        if len(pts) > 3:
            pts = pts[ConvexHull(pts).vertices]
    except Exception:
        pass
    if len(pts) > max_points:
        idx = np.random.default_rng(0).choice(len(pts), max_points, replace=False)
        pts = pts[idx]
    return float(pdist(pts).max())


def isolate_attached_flame(mask):
    """Keep only the LARGEST connected component of `mask` (by cell
    count) -- empirically this is always the real, continuously-fed,
    nozzle-attached flame (tens of thousands of cells), while spurious
    disconnected blobs (a handful to a few thousand cells) are numerical
    noise or transients. In particular, this case briefly injects a hot
    ignition pulse at the nozzle at simulation start
    (prob.start_time_injProd / duration_injProd in input.3d-regt); as
    that slug is swept downstream and out through the outflow it forms a
    T>=800K region fully DISCONNECTED from the attached flame, which
    would otherwise corrupt the "farthest downstream" length measurement.
    (Picking the component with smallest y instead of largest size was
    tried and rejected: it gets fooled by single stray hot cells sitting
    right at/near the nozzle.)"""
    labeled, n = label(mask)
    if n <= 1:
        return mask
    counts = np.bincount(labeled.ravel())
    counts[0] = 0  # background is never the flame
    return labeled == np.argmax(counts)


def flame_length_width(T, x, y, z, volfrac=None):
    """T: 3D array [Nx,Ny,Nz] on the uniform grid defined by the 1D
    coordinate arrays x,y,z (meters). Returns (length, width) in meters,
    both NaN if no cell reaches T_FLAME."""
    mask = T >= T_FLAME
    if volfrac is not None:
        mask &= volfrac > 1e-3  # exclude EB-covered/solid nozzle cells

    if not np.any(mask):
        return np.nan, np.nan

    mask = isolate_attached_flame(mask)

    y_reach = y[np.any(mask, axis=(0, 2))]
    length = float(y_reach.max() - JET_Y)

    footprint = np.any(mask, axis=1)  # (Nx,Nz): envelope projected along y
    ix, iz = np.nonzero(footprint)
    width = enclosing_diameter(x[ix], z[iz])

    return length, width


def uniform_grid_and_temperature(ds, level):
    dims = (ds.domain_dimensions * ds.refine_by**level).astype(int)
    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)
    T = cg["temp"].v
    has_volfrac = ("boxlib", "volFrac") in ds.field_list
    vfrac = cg["volFrac"].v if has_volfrac else None

    le = ds.domain_left_edge.v
    re_ = ds.domain_right_edge.v
    dx = (re_ - le) / dims
    x = le[0] + (np.arange(dims[0]) + 0.5) * dx[0]
    y = le[1] + (np.arange(dims[1]) + 0.5) * dx[1]
    z = le[2] + (np.arange(dims[2]) + 0.5) * dx[2]
    return T, vfrac, x, y, z


def plot_midplane(T_mean, x, y, z, length_avg, width_avg, plotfiles, out_path):
    k_mid = int(np.argmin(np.abs(z - 0.5 * (z[0] + z[-1]))))
    slice_xy = T_mean[:, :, k_mid]  # (Nx, Ny)

    fig, ax = plt.subplots(figsize=(11, 5))
    im = ax.imshow(
        slice_xy, origin="lower", aspect="auto", cmap="inferno",
        extent=[y[0], y[-1], x[0], x[-1]],
    )
    ax.contour(y, x, slice_xy, levels=[T_FLAME], colors="cyan", linewidths=1.5)
    ax.axvline(JET_Y, color="white", linestyle=":", linewidth=1, label="Nozzle (y)")
    ax.axvline(JET_Y + length_avg, color="lime", linestyle="--", linewidth=1.5,
               label=f"Length = {length_avg:.3f} m")
    ax.set_xlabel("y (downstream) [m]")
    ax.set_ylabel("x (jet direction) [m]")
    ax.set_title(f"Time-averaged T=800K flame envelope, z-midplane\n"
                 f"{plotfiles[0]} - {plotfiles[-1]}  |  Length={length_avg:.3f} m, "
                 f"Width={width_avg:.3f} m")
    fig.colorbar(im, ax=ax, label="Temperature [K]")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    print(f"Saved averaged-flame midplane plot to {out_path}")


def main():
    first_plt = sys.argv[1] if len(sys.argv) > 1 else FIRST_PLT
    last_plt = sys.argv[2] if len(sys.argv) > 2 else LAST_PLT

    plotfiles = get_plotfiles_in_range(first_plt, last_plt)
    if not plotfiles:
        print(f"No plt directories found between {first_plt} and {last_plt}.")
        return
    print(f"Found {len(plotfiles)} plotfiles: {plotfiles[0]} ... {plotfiles[-1]}")

    os.makedirs(OUT_DIR, exist_ok=True)

    # Fix the resampling level once (from the first file) so every snapshot's
    # covering_grid has the SAME shape and can be accumulated into a mean field.
    level = min(GRID_LEVEL, yt.load(plotfiles[0]).max_level)
    print(f"Resampling every snapshot onto AMR level {level} for the shape analysis.")

    times, lengths, widths = [], [], []
    T_sum = Vf_sum = x_ref = y_ref = z_ref = None
    n_avg = 0

    for pf in plotfiles:
        print(f"Processing {pf} ...")
        try:
            ds = yt.load(pf)
            T, vfrac, x, y, z = uniform_grid_and_temperature(ds, level)
        except Exception as e:
            print(f"  [WARN] skipping {pf}: {e}")
            continue

        t = ds.current_time.to_value("s")
        length, width = flame_length_width(T, x, y, z, vfrac)
        times.append(t)
        lengths.append(length)
        widths.append(width)
        print(f"  t = {t:.4f} s | length = {length:.4f} m ({length / D_NOZZLE:.2f} D) "
              f"| width = {width:.4f} m ({width / D_NOZZLE:.2f} D)")

        if T_sum is None:
            T_sum = np.zeros_like(T)
            Vf_sum = np.zeros_like(T) if vfrac is not None else None
            x_ref, y_ref, z_ref = x, y, z
        T_sum += T
        if Vf_sum is not None:
            Vf_sum += vfrac
        n_avg += 1

    if not times:
        print("No plotfile could be processed.")
        return

    times = np.array(times)
    lengths = np.array(lengths)
    widths = np.array(widths)

    # ---- averaged-field result (paper-equivalent metric) ----
    T_mean = T_sum / n_avg
    Vf_mean = (Vf_sum / n_avg) if Vf_sum is not None else None
    length_avg, width_avg = flame_length_width(T_mean, x_ref, y_ref, z_ref, Vf_mean)

    print("\n--- TIME-AVERAGED FLAME SHAPE (Bradley et al. 2016 / Mohit et al. 2026 method) ---")
    print(f"Averaged over {n_avg} snapshots, t = {times.min():.4f} - {times.max():.4f} s")
    print(f"Length = {length_avg:.4f} m ({length_avg / D_NOZZLE:.2f} D)")
    print(f"Width  = {width_avg:.4f} m ({width_avg / D_NOZZLE:.2f} D)")

    # ---- save time series ----
    df = pd.DataFrame({"time_s": times, "length_m": lengths, "width_m": widths})
    csv_path = os.path.join(OUT_DIR, "flame_shape_history.csv")
    df.to_csv(csv_path, index=False)
    print(f"\nSaved time history to {csv_path}")

    summary_path = os.path.join(OUT_DIR, "flame_shape_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Flame shape summary (T=800K envelope, Mohit et al. 2026 method)\n")
        f.write(f"Plotfiles: {plotfiles[0]} .. {plotfiles[-1]} ({n_avg} snapshots)\n")
        f.write(f"Time range: {times.min():.4f} - {times.max():.4f} s\n")
        f.write(f"AMR level used: {level}\n")
        f.write(f"Nozzle diameter D = {D_NOZZLE:.4f} m, length origin y_jet = {JET_Y:.4f} m\n")
        f.write(f"Time-averaged length = {length_avg:.4f} m ({length_avg / D_NOZZLE:.2f} D)\n")
        f.write(f"Time-averaged width  = {width_avg:.4f} m ({width_avg / D_NOZZLE:.2f} D)\n")
        f.write(f"Instantaneous length: mean={lengths.mean():.4f} m, std={lengths.std():.4f} m\n")
        f.write(f"Instantaneous width:  mean={widths.mean():.4f} m, std={widths.std():.4f} m\n")
    print(f"Saved summary to {summary_path}")

    # ---- time-history plots ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax1.plot(times, widths, color="tab:blue", marker="o", linestyle="-",
             linewidth=1.5, label="Instantaneous width")
    ax1.axhline(width_avg, color="tab:blue", linestyle="--", alpha=0.7,
                label=f"Time-averaged width ({width_avg:.3f} m)")
    ax1.set_ylabel("Flame width [m]", fontsize=12)
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="best")

    ax2.plot(times, lengths, color="tab:red", marker="o", linestyle="-",
             linewidth=1.5, label="Instantaneous length")
    ax2.axhline(length_avg, color="tab:red", linestyle="--", alpha=0.7,
                label=f"Time-averaged length ({length_avg:.3f} m)")
    ax2.set_xlabel("Time [s]", fontsize=12)
    ax2.set_ylabel("Flame length [m]", fontsize=12)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="best")

    fig.suptitle(f"Flame shape history (T = {T_FLAME:.0f} K envelope), "
                 f"{plotfiles[0]}-{plotfiles[-1]}")
    fig.tight_layout()
    plot_path = os.path.join(OUT_DIR, "flame_shape_history.png")
    fig.savefig(plot_path, dpi=300)
    print(f"Saved width/length history plot to {plot_path}")

    # ---- averaged-envelope visualization (bonus, mirrors paper's Fig. 5) ----
    midplane_path = os.path.join(OUT_DIR, "flame_shape_averaged_midplane.png")
    plot_midplane(T_mean, x_ref, y_ref, z_ref, length_avg, width_avg, plotfiles, midplane_path)


if __name__ == "__main__":
    main()
