#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import numpy as np
import os
import sys 
import glob

# Global settings
filename = "amrex"

def get_safe_extent(x, y, default=[0, 1, 0, 1]):
    if np.all(np.isnan(x)) or np.all(np.isnan(y)):
        return default
    return [np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)]

def plot_rho(plt_path, output_prefix=""):
    """
    Loads a single plt file, plots data with AMReX mesh lines, and saves to a PNG.
    """
    # Load plotfile
    ds = yt.load(plt_path)
    
    # Define field
    #field_name = ('boxlib', 'MANI_PROG')
    field_name = ('boxlib', 'I_R(PROG)')
    #field_name = ('boxlib', 'Y(ZMIX)')
    #field_name = ('boxlib', 'temp')
    zmix_field = ('boxlib', 'Y(ZMIX)')
    
    if field_name not in ds.field_list:
        print(f"  [ERROR] Field '{field_name}' not found in {plt_path}. Skipping.")
        return

    # Extract data as a fixed-resolution buffer (FRB)
    slc = yt.SlicePlot(ds, "z", field_name)
    frb = slc.frb
    
    # Get coordinates and shift to start at 0,0 based on domain edge
    x_coords = frb["x"].v - ds.domain_left_edge[0].v
    y_coords = frb["y"].v - ds.domain_left_edge[1].v
    data = frb[field_name].v
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot the scalar field
    img = ax.imshow(
        data,
        extent=get_safe_extent(x_coords, y_coords),
        origin="lower",
        cmap="RdBu_r",
        aspect="equal",
    )
    
    # --- STOICHIOMETRIC ISOLINE OVERLAY ---
    if zmix_field in ds.field_list:
        zmix_data = frb[zmix_field].v
        # Plot the 0.055 contour line
        iso_contour = ax.contour(
            x_coords, y_coords, zmix_data,
            levels=[0.055],
            colors='black',       # High contrast line color
            linewidths=1.5,       # Clean visible thickness
            linestyles='solid'
        )
        # Optional: Add a text label directly on the isoline
        ax.clabel(iso_contour, inline=True, fmt='Z=0.055', fontsize=9, colors='black')
    else:
        print(f"  [WARNING] Iso-contour field '{zmix_field}' not found in {plt_path}. Skipping line.")
    
    # --- AMReX MESH OVERLAY ---
    # Define colors for each level: Level 0 -> black, Level 1 -> red, Level 2 -> green, etc.
    level_colors = ['black', 'red', 'green', 'blue', 'darkorange', 'purple', 'cyan']
    
    # Keep track of which levels have been plotted to display in the legend
    levels_plotted = set()
    
    # Loop through all grids in the AMR hierarchy
    for grid in ds.index.grids:
        # Get the AMR level integer (0, 1, 2, etc.)
        lvl = grid.Level
        levels_plotted.add(lvl)
        
        # Pick the color. If level exceeds the predefined list, fallback safely using modulo
        color = level_colors[lvl % len(level_colors)]
        
        # Calculate edges relative to the domain origin (matching the imshow shift)
        lower_left = grid.LeftEdge.v - ds.domain_left_edge.v
        upper_right = grid.RightEdge.v - ds.domain_left_edge.v
        
        width = upper_right[0] - lower_left[0]
        height = upper_right[1] - lower_left[1]
        
        # Create a rectangle for the grid boundary with level-specific color
        rect = patches.Rectangle(
            (lower_left[0], lower_left[1]), 
            width, 
            height,
            linewidth=0.8,       # Slightly increased width to make colored lines pop
            edgecolor=color,     # Dynamic level color
            facecolor='none', 
            alpha=0.6            # Transparency for better visibility of data
        )
        ax.add_patch(rect)
    # ---------------------------

    # --- ADD LEGEND ---
    # Create custom legend elements based on the levels we actually encountered
    legend_elements = [
        Line2D([0], [0], color=level_colors[lvl % len(level_colors)], lw=2, label=f'Level {lvl}')
        for lvl in sorted(levels_plotted)
    ]
    
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', title='AMR Levels')
    # ------------------

    # Formatting
    plt.colorbar(img, ax=ax, label=field_name[1])
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    
    plot_name = os.path.basename(plt_path)
    
    # Save
    save_path = f"{output_prefix}.png"
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot AMReX data with mesh lines.')
    parser.add_argument("base_dir", help="Path to directory containing pltXXXXX folders")
    parser.add_argument("-o", "--output", default=filename, help="Output filename prefix")
    args = parser.parse_args()

    plot_pattern = os.path.join(args.base_dir, "plt*")
    plot_dirs = sorted([f for f in glob.glob(plot_pattern) if os.path.isdir(f)])

    if not plot_dirs:
        print(f"Error: No 'plt*' directories found in {args.base_dir}")
        sys.exit(1)

    print(f"--- Found {len(plot_dirs)} plotfiles to process ---")

    for plt_path in plot_dirs:
        plot_name = os.path.basename(plt_path)
        output_save_path = os.path.join(args.base_dir, f"{args.output}_{plot_name}.png")
        output_prefix_for_func = os.path.join(args.base_dir, f"{args.output}_{plot_name}")

        if os.path.exists(output_save_path):
            print(f"Skipping {plot_name} (image already exists).")
            continue

        try:
            print(f"Processing: {plot_name}")
            plot_rho(plt_path, output_prefix_for_func)
        except Exception as e:
            print(f"  [ERROR] Failed to process {plot_name}: {e}")
    
    print("\n--- Processing complete ---")