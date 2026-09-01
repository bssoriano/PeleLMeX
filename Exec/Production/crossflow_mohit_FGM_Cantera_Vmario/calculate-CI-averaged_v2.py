#!/usr/bin/env python3
"""
Combustion inefficiency (CI = 1 - eta) history for the crossflow FGM case,
following Mohit et al. (2026, Flow Turbul. Combust. 116:6), Sec. 2.1, Eq.
(2)-(3):

    eta = mdot_C,as_CO2,out / mdot_C,in     (fraction of the fuel jet's
                                              carbon that leaves as CO2)
    CI  = 1 - eta

mdot_C,in = rho_j * A_j * U_j * Y_C is fixed by the (constant) jet BC.
mdot_C,out is measured here by the mass-conservation COMPLEMENT of Eq. (3):
instead of integrating the CO2 flux directly, this integrates the UNBURNED
carbon flux (CH4 + CO -- the only other carbon-bearing species with
non-negligible mass in this manifold table; CH2O was checked directly
against a plt file and found ~1e4x smaller) leaving through the domain's
one outflow plane. By carbon conservation this equals the paper's CI
exactly (same number, complementary route), and only needs species that
are already in the plt files.

Fixes applied after comparing the original version of this script against
the paper (2026-08-09 chat):
  1. Jet density now uses the ACTUAL prob.jet_T from input.3d-regt (parsed
     at runtime, not hardcoded) instead of a stale "298 K" guess -- that
     value was actually prob.T_mean (the ambient/crossflow temperature),
     not the jet's. This alone changes mdot_C,in by ~2.3x.
  2. Nozzle area kept as a FULL circle (pi*r^2): confirmed correct for
     this case's protruding EB nozzle (fuel exits the pipe's flat tip,
     which is not bisected by any domain wall) -- unlike calculate-CI.py's
     0.5*pi*r^2, which looks left over from an older flush-nozzle variant
     and is NOT fixed here since that's a different file.
  3. Outlet sampling plane moved from an arbitrary 90%/95%-of-domain
     fraction to the LAST valid interior cell layer next to the true
     outflow boundary (y = domain_right_edge - dy/2). Checked empirically:
     CH4/CO/CO2 flux differs by 20-35% between the old 90% and 95% planes,
     so this choice matters and should match the paper's "at the exit
     plane" definition as closely as a cell-centered AMR field allows.
Also added: flow-through-time / paper-equivalent averaging-window context
and a first-half-vs-second-half stationarity ("plateau") check, mirroring
flame_shape_history.py -- CI is only meaningful to average over a
statistically steady window, and this case's own prior output
(historial_combustion_separado.png) shows CI still ramping 0%->25% with no
plateau over the only data available so far.

Usage:
    python3 calculate-CI-averaged.py [first_plt] [last_plt]
Both arguments are optional; if omitted, FIRST_PLT/LAST_PLT below are used.
"""
import matplotlib
matplotlib.use("Agg")

import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yt

yt.funcs.mylog.setLevel(50)

# ---------------------------------------------------------------------
# CONFIGURATION - edit these, or override first/last plt on the command line
# ---------------------------------------------------------------------
FIRST_PLT = "plt00000"
LAST_PLT =  "plt133981"

INPUT_FILE = "input.3d-regt"   # jet/crossflow BC parameters are read from here
OUT_DIR = "combustion_efficiency_results"
# ---------------------------------------------------------------------


def read_prob_param(input_file, key, cast=float):
    """Pull a `key = value` line (e.g. 'prob.jet_T') straight out of the
    AMReX ParmParse input file, so this script can't silently drift out of
    sync with the actual simulation setup (as the hardcoded T_jet did)."""
    value = None
    with open(input_file) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            if k.strip() == key:
                value = cast(v.strip().split()[0])  # last match wins, like ParmParse
    if value is None:
        raise KeyError(f"'{key}' not found in {input_file}")
    return value


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


def main():
    first_plt = sys.argv[1] if len(sys.argv) > 1 else FIRST_PLT
    last_plt = sys.argv[2] if len(sys.argv) > 2 else LAST_PLT

    plotfiles = get_plotfiles_in_range(first_plt, last_plt)
    if not plotfiles:
        print(f"No plt directories found between {first_plt} and {last_plt}.")
        return
    print(f"Found {len(plotfiles)} plotfiles: {plotfiles[0]} ... {plotfiles[-1]}")
    os.makedirs(OUT_DIR, exist_ok=True)

    # -------------------------------------------------------------
    # 1. INLET CARBON FLOW RATE (Eq. 2 denominator), from the ACTUAL jet
    #    BC parameters in input.3d-regt -- not hardcoded guesses (fix #1).
    # -------------------------------------------------------------
    W_C, W_CH4, W_CO = 12.011, 16.042, 28.010
    Z_C_CH4 = W_C / W_CH4
    Z_C_CO = W_C / W_CO
    R_univ = 8.31446  # J/(mol K)

    P_mean = read_prob_param(INPUT_FILE, "prob.P_mean")     # Pa
    T_jet = read_prob_param(INPUT_FILE, "prob.jet_T")        # K
    v_jet = read_prob_param(INPUT_FILE, "prob.jet_vel")      # m/s
    r_jet = read_prob_param(INPUT_FILE, "prob.jet_radius")   # m
    Uc = read_prob_param(INPUT_FILE, "prob.meanFlowMag")     # m/s

    rho_jet = (P_mean * (W_CH4 / 1000.0)) / (R_univ * T_jet)
    A_jet = np.pi * r_jet**2   # full circle (fix #2, see module docstring)
    flujo_CH4_in = rho_jet * A_jet * v_jet
    flujo_C_in = flujo_CH4_in * Z_C_CH4

    print("\n--- INLET CARBON FLOW (from input.3d-regt, Eq. 2 of Mohit et al. 2026) ---")
    print(f"P_mean={P_mean:.1f} Pa, T_jet={T_jet:.1f} K, v_jet={v_jet:.3f} m/s, r_jet={r_jet:.4f} m")
    print(f"rho_jet = {rho_jet:.4f} kg/m3, A_jet = {A_jet:.6e} m2")
    print(f"mdot_CH4_in = {flujo_CH4_in:.6e} kg/s")
    print(f"mdot_C_in   = {flujo_C_in:.6e} kg/s")

    # -------------------------------------------------------------
    # 2. OUTLET CARBON FLOW RATE (unburned CH4+CO complement of Eq. 3),
    #    sampled at the LAST interior cell layer next to the true outflow
    #    boundary, not an arbitrary 90%/95% offset (fix #3).
    # -------------------------------------------------------------
    times, ineficiencias, eficiencias = [], [], []
    domain_left = domain_right = None

    print("\n--- PROCESSING TIME HISTORY (PeleLMeX) ---")
    for pf in plotfiles:
        print(f"Loading {pf}...")
        try:
            ds = yt.load(pf)
        except Exception as e:
            print(f"  [ERROR] could not load {pf}: {e}")
            continue

        if domain_left is None:
            domain_left = ds.domain_left_edge.v
            domain_right = ds.domain_right_edge.v

        t_sim = ds.current_time.to_value("s")
        dy = ds.index.get_smallest_dx().v
        y_outlet = ds.domain_right_edge[1].v - dy / 2.0

        try:
            slc = ds.slice("y", y_outlet)
            rho = slc[("boxlib", "density")].v
            vy = slc[("boxlib", "y_velocity")].v
            dx = slc[("index", "dx")].v
            dz = slc[("index", "dz")].v

            field_ch4 = ("boxlib", "MANI_Y-CH4")
            field_co = ("boxlib", "MANI_Y-CO")
            flujo_CH4_out = np.sum(rho * vy * slc[field_ch4].v * dx * dz) if field_ch4 in ds.field_list else 0.0
            flujo_CO_out = np.sum(rho * vy * slc[field_co].v * dx * dz) if field_co in ds.field_list else 0.0
        except Exception as e:
            print(f"  [ERROR] integration failed for {pf}: {e}")
            continue

        flujo_C_out = flujo_CH4_out * Z_C_CH4 + flujo_CO_out * Z_C_CO
        ineficiencia = (flujo_C_out / flujo_C_in) * 100.0
        eficiencia = 100.0 - ineficiencia

        times.append(t_sim)
        ineficiencias.append(ineficiencia)
        eficiencias.append(eficiencia)
        print(f"  t = {t_sim:.4f} s | CI = {ineficiencia:.4f} % | eta = {eficiencia:.4f} %")

    if not times:
        print("No plotfile could be processed.")
        return

    times = np.array(times)
    ineficiencias = np.array(ineficiencias)
    eficiencias = np.array(eficiencias)

    avg_ineficiencia = np.mean(ineficiencias)
    avg_eficiencia = np.mean(eficiencias)

    # -------------------------------------------------------------
    # 3. STATIONARITY / "HAS-IT-PLATEAUED" CHECK. CI is only meaningful
    #    to average over a statistically steady window; the paper fixes
    #    tau=8s ~ 630 D/Uc AFTER reaching steady state. This case's
    #    analogous window is computed below for reference.
    # -------------------------------------------------------------
    D_nozzle = 2.0 * r_jet
    y_len = domain_right[1] - domain_left[1]
    flow_through_time = y_len / Uc
    paper_equiv_window = 630.0 * D_nozzle / Uc

    half = len(times) // 2
    if half >= 2:
        mean1, mean2 = ineficiencias[:half].mean(), ineficiencias[half:].mean()
        std1, std2 = ineficiencias[:half].std(), ineficiencias[half:].std()
        plateaued = abs(mean2 - mean1) < (std1 + std2)
    else:
        mean1 = mean2 = std1 = std2 = float("nan")
        plateaued = False

    print("\n--- STATIONARITY CHECK ---")
    print(f"Simulated window: t = {times.min():.4f} - {times.max():.4f} s "
          f"(duration {times.max() - times.min():.4f} s)")
    print(f"Domain flow-through time (downstream length / Uc) = {flow_through_time:.4f} s")
    print(f"Paper-equivalent averaging window (630 D/Uc, this case's D & Uc) = {paper_equiv_window:.4f} s")
    print(f"First-half CI: {mean1:.3f}% +/- {std1:.3f}  |  Second-half CI: {mean2:.3f}% +/- {std2:.3f}")
    if plateaued:
        print("=> First/second-half means agree within their fluctuation level: window looks plateaued.")
    else:
        print("=> [WARNING] first/second-half means differ by more than their fluctuation level: this "
              "window is still drifting (likely a startup transient), NOT a converged average yet. "
              "Re-run with a later start plt once the CI-vs-time curve visibly flattens.")

    print("\n--- TIME-AVERAGED RESULT ---")
    print(f"Combustion Inefficiency (CI): {avg_ineficiencia:.6f} %")
    print(f"Combustion Efficiency (eta):  {avg_eficiencia:.6f} %")

    # -------------------------------------------------------------
    # 4. SAVE + PLOT
    # -------------------------------------------------------------
    df = pd.DataFrame({"time_s": times, "CI_percent": ineficiencias, "eta_percent": eficiencias})
    csv_path = os.path.join(OUT_DIR, "combustion_efficiency_history.csv")
    df.to_csv(csv_path, index=False)

    summary_path = os.path.join(OUT_DIR, "combustion_efficiency_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Combustion inefficiency summary (Mohit et al. 2026, Eq. 2-3 equivalent)\n")
        f.write(f"Plotfiles: {plotfiles[0]} .. {plotfiles[-1]} ({len(times)} snapshots)\n")
        f.write(f"T_jet={T_jet}K, P_mean={P_mean}Pa, v_jet={v_jet}m/s, r_jet={r_jet}m (from {INPUT_FILE})\n")
        f.write(f"mdot_C_in = {flujo_C_in:.6e} kg/s\n")
        f.write(f"Time-averaged CI  = {avg_ineficiencia:.4f} %\n")
        f.write(f"Time-averaged eta = {avg_eficiencia:.4f} %\n")
        f.write(f"Flow-through time = {flow_through_time:.4f} s; paper-equivalent window = {paper_equiv_window:.4f} s\n")
        f.write(f"Plateaued: {plateaued} (first-half {mean1:.3f}%, second-half {mean2:.3f}%)\n")
    print(f"\nSaved {csv_path}")
    print(f"Saved {summary_path}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    color_inef = "tab:red"
    ax1.set_ylabel("Combustion Inefficiency, CI (%)", color=color_inef, fontsize=12)
    ax1.plot(times, ineficiencias, color=color_inef, marker="o", linestyle="-", linewidth=2, label="CI (instantaneous)")
    ax1.tick_params(axis="y", labelcolor=color_inef)
    ax1.axhline(avg_ineficiencia, color=color_inef, linestyle="--", alpha=0.6,
                label=f"Time-averaged CI ({avg_ineficiencia:.3f}%)")
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="best")
    flag = "plateaued" if plateaued else "NOT YET PLATEAUED -- treat average as preliminary"
    ax1.set_title(f"Combustion inefficiency history ({flag})", fontsize=12)

    color_ef = "tab:blue"
    ax2.set_xlabel("Time (s)", fontsize=12)
    ax2.set_ylabel("Combustion Efficiency (%)", color=color_ef, fontsize=12)
    ax2.plot(times, eficiencias, color=color_ef, marker="x", linestyle="-", linewidth=1.5, label="Efficiency")
    ax2.tick_params(axis="y", labelcolor=color_ef)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="best")

    fig.tight_layout()
    plot_path = os.path.join(OUT_DIR, "combustion_efficiency_history.png")
    fig.savefig(plot_path, dpi=300)
    print(f"Saved plot to {plot_path}")


if __name__ == "__main__":
    main()
