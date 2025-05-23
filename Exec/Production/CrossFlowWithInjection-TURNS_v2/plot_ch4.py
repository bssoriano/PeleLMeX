#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np

def plot_CH4_massfrac(plt_path, output_prefix="CH4"):
    # Load plotfile
    ds = yt.load(plt_path)
    
    # Extract data as a fixed-resolution buffer (FRB)
    slc = yt.SlicePlot(ds, "z", ("boxlib", "Y(CH4)"))
    frb = slc.frb
    
    # Get coordinates (shift to start at 0,0)
    x = frb["x"].v - ds.domain_left_edge[0].v  # Shift to origin
    y = frb["y"].v - ds.domain_left_edge[1].v  # Shift to origin
    CH4 = frb[("boxlib", "Y(CH4)")].v
    
    # Create figure
    plt.figure(figsize=(8, 8))  # Square figure for 1:1 aspect
    img = plt.imshow(
        CH4,
        extent=[x.min(), x.max(), y.min(), y.max()],
        origin="lower",  # Origin at bottom-left
        cmap="RdBu_r",
        aspect="equal",  # Enforce 1:1 aspect ratio
    )
    plt.colorbar(label="CH4 Mass Fraction")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("CH4 Mass Fraction")
    
    # Save
    plt.savefig(f"{output_prefix}.png", bbox_inches="tight", dpi=300)
    print(f"Saved plot to {output_prefix}.png")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Plot Y(CH4) with origin at (0,0) and 1:1 aspect')
    parser.add_argument("plt_dir", help="Path to pltXXXXX directory")
    parser.add_argument("-o", "--output", default="CH4", help="Output filename prefix")
    args = parser.parse_args()
    plot_CH4_massfrac(args.plt_dir, args.output)