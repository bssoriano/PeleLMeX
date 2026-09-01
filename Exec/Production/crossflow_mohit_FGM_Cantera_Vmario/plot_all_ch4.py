#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
import glob
import argparse

filename = "zmix_plot"

def plot_temperature(plt_path, output_prefix=""):
    """
    Loads a 2D PeleLMex plotfile, extracts temperature and velocity 
    in the XY plane, and saves to a PNG using matplotlib with streamlines.
    """
    # Load plotfile
    ds = yt.load(plt_path)
    
    # PeleLMex fields
    #field_temp = ('boxlib', 'temp')
    field_temp = ('boxlib', 'Y(ZMIX)')
    #field_temp = ('boxlib', 'Y(PROG)')
    field_u = ('boxlib', 'x_velocity')
    field_v = ('boxlib', 'y_velocity')
    
    # Check if the temperature field actually exists in the plotfile
    if field_temp not in ds.field_list:
        print(f"  [ERROR] Field '{field_temp}' not found in {plt_path}. Skipping.")
        return

    # Slice along "z" and extract all three fields at once
    slc = yt.SlicePlot(ds, "z", [field_temp, field_u, field_v])
    frb = slc.frb
    
    # Extract the bare numpy arrays for Temp, U, and V
    data_temp = frb[field_temp].v
    data_u = frb[field_u].v
    data_v = frb[field_v].v

    # Calculate domain bounds and shift to origin (0,0)
    x_min = ds.domain_left_edge[0].v
    x_max = ds.domain_right_edge[0].v
    y_min = ds.domain_left_edge[1].v
    y_max = ds.domain_right_edge[1].v
    
    # extent = [left, right, bottom, top]
    extent = [x_min, x_max, y_min, y_max]    
    
    # --- Generate Coordinate Grids for Streamlines ---
    # Streamplot needs X and Y arrays that match the shape of our data arrays.
    ny, nx = data_temp.shape
    x_coords = np.linspace(x_min, x_max, nx)
    y_coords = np.linspace(y_min, y_max, ny)
    X, Y = np.meshgrid(x_coords, y_coords)
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot the temperature background
    img = ax.imshow(
        data_temp,
        extent=extent,
        origin="lower",  
        cmap="RdBu_r",  # "inferno" or "hot" are standard for temperature 
        aspect="equal",  
    )
    
    plt.colorbar(img, ax=ax, label="Mixture Fraction")
    
    # --- Add the Streamlines ---
    # We plot the U and V arrays over our coordinate meshgrid.
    # Adjust 'density' to control how tightly packed the streamlines are.
    # Adjust 'linewidth' or 'color' as needed.
    ax.streamplot(X, Y, data_u, data_v, color='black', linewidth=0.8, density=1.5)
      
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    
    # Add a title with the plotfile name
    plot_name = os.path.basename(plt_path)
    ax.set_title(f"Temperature and Flow Field - {plot_name}")
    
    # Save
    save_path = f"{output_prefix}.png"
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    
    # CRITICAL: Close the figure to free memory
    plt.close(fig)

# --- Main execution block ---
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plot temperature and streamlines for all timesteps in a directory.')
    parser.add_argument("base_dir", help="Path to the directory containing pltXXXXX folders")
    parser.add_argument("-o", "--output", default=filename, help="Output filename prefix (e.g., 'temp_plot')")
    args = parser.parse_args()

    # Find all directories matching "plt*" in the base directory
    plot_pattern = os.path.join(args.base_dir, "plt*")
    plot_dirs = sorted([f for f in glob.glob(plot_pattern) if os.path.isdir(f)])

    if not plot_dirs:
        print(f"Error: No 'plt*' directories found in {os.path.abspath(args.base_dir)}")
        sys.exit(1)

    print(f"--- Found {len(plot_dirs)} plotfiles to process ---")

    # Loop over each plot directory
    for plt_path in plot_dirs:
        plot_name = os.path.basename(plt_path)
        output_filename = f"{args.output}_{plot_name}.png"
        output_save_path = os.path.join(args.base_dir, output_filename)
        output_prefix_for_func = os.path.join(args.base_dir, f"{args.output}_{plot_name}")

        if os.path.exists(output_save_path):
            print(f"Skipping {plot_name} (image already exists).")
            continue

        try:
            print(f"Processing: {plot_name} -> {output_filename}")
            plot_temperature(plt_path, output_prefix_for_func)
        except Exception as e:
            print(f"  [ERROR] Failed to process {plot_name}: {e}")
    
    print("\n--- Processing complete ---")