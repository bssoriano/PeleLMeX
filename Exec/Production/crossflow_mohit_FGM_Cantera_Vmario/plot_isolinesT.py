#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
import glob

# --- PARAMETERS ---
D = 0.0762             # Diameter in meters
iso_temp = 800         # Temperature for the isoline [K]
filename = f"temp_isoline_{iso_temp}K"

def get_safe_extent(x, y, default=[0, 1, 0, 1]):
    if np.all(np.isnan(x)) or np.all(np.isnan(y)):
        return default
    return [np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)]

def plot_rho(plt_path, output_prefix=""):
    """
    Loads a single plt file, plots a temperature isoline on the yz plane, and saves to a PNG.
    """
    # Load plotfile
    ds = yt.load(plt_path)
    
    # --- CHANGED: Field to plot is now temperature ---
    field_name = ('boxlib', 'temp')
    
    if field_name not in ds.field_list:
        print(f"  [ERROR] Field '{field_name}' not found in {plt_path}. Skipping.")
        return

    # Extract data on the YZ plane at x = 0.381
    center = [0.381, ds.domain_center[1].v, ds.domain_center[2].v]
    slc = yt.SlicePlot(ds, "x", field_name, center=center)
    frb = slc.frb
    
    # Non-dimensionalize the coordinates
    y_coords = (frb["y"].v - ds.domain_left_edge[1].v) / D
    z_coords = frb["z"].v / D 
    
    # Transpose data so Z is on horizontal, Y is on vertical
    data = frb[field_name].v.T
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(8, 8))  # Square figure
    
    # Get the min/max of our new non-dimensional coordinates
    z_min, z_max = np.nanmin(z_coords), np.nanmax(z_coords)
    y_min, y_max = np.nanmin(y_coords), np.nanmax(y_coords)
    
    # --- CHANGED: Replaced imshow with contour for the isoline ---
    # Generate 1D arrays and a meshgrid to match the data shape for contouring
    z_1d = np.linspace(z_min, z_max, data.shape[1])
    y_1d = np.linspace(y_min, y_max, data.shape[0])
    ZZ, YY = np.meshgrid(z_1d, y_1d)
    
    # Plot the isoline if the temperature exists in this slice
    if np.nanmin(data) <= iso_temp <= np.nanmax(data):
        CS = ax.contour(ZZ, YY, data, levels=[iso_temp], colors='red', linewidths=2)
    #    ax.clabel(CS, inline=True, fontsize=12, fmt='%1.0f K')
    else:
        print(f"  [WARNING] {iso_temp} K not found in {os.path.basename(plt_path)}.")
    
    # Add the circle: Center at z/D = 0, y/D = 10. Radius = 0.5
    circle = plt.Circle((0.0, 10.0), radius=0.5, color='black', fill=False, linestyle='-', linewidth=2, label='D=1')
    ax.add_patch(circle)
    
    # Update axis labels and enforce equal aspect ratio (since imshow is gone)
    ax.set_xlabel("z / D")
    ax.set_ylabel("y / D")
    ax.set_aspect('equal', adjustable='box')

    ax.set_xlim(-2,2)
    ax.set_ylim(8,16)
    
    plot_name = os.path.basename(plt_path)
    ax.set_title(f"{plot_name} | Temp Isoline ({iso_temp} K) | Plane yz at x = 0.381 m")
    
    # Save
    save_path = f"{output_prefix}.png"
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

# --- Main execution block ---
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot field for all timesteps in a directory.')
    parser.add_argument("base_dir", help="Path to the directory containing pltXXXXX folders")
    parser.add_argument("-o", "--output", default=filename, help="Output filename prefix")
    args = parser.parse_args()

    plot_pattern = os.path.join(args.base_dir, "plt*")
    plot_dirs = sorted([f for f in glob.glob(plot_pattern) if os.path.isdir(f)])

    if not plot_dirs:
        print(f"Error: No 'plt*' directories found in {os.path.abspath(args.base_dir)}")
        sys.exit(1)

    print(f"--- Found {len(plot_dirs)} plotfiles to process ---")

    # Create the "figures" directory
    fig_dir = "figures"
    os.makedirs(fig_dir, exist_ok=True)
    print(f"Figures will be saved to: ./{fig_dir}/")

    for plt_path in plot_dirs:
        plot_name = os.path.basename(plt_path)
        output_filename = f"{args.output}_{plot_name}.png"
        
        # Route the saves to the 'figures' folder
        output_save_path = os.path.join(fig_dir, output_filename)
        output_prefix_for_func = os.path.join(fig_dir, f"{args.output}_{plot_name}")

        if os.path.exists(output_save_path):
            print(f"Skipping {plot_name} (image already exists).")
            continue

        try:
            print(f"Processing: {plot_name} -> {output_filename}")
            plot_rho(plt_path, output_prefix_for_func)
        except Exception as e:
            print(f"  [ERROR] Failed to process {plot_name}: {e}")
    
    print("\n--- Processing complete ---")