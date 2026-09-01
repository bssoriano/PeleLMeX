#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
import glob

#filename = "mag_vel"

filename = "temp"

def get_safe_extent(x, y, default=[0, 1, 0, 1]):
    # Check if arrays are empty or all-NaN
    if np.all(np.isnan(x)) or np.all(np.isnan(y)):
        return default

    return [
        np.nanmin(x), 
        np.nanmax(x), 
        np.nanmin(y), 
        np.nanmax(y)
    ]

def plot_rho(plt_path, output_prefix=""):
    """
    Loads a single plt file, plots density, and saves to a PNG.
    """
    # Load plotfile
    ds = yt.load(plt_path)
    
    # The standard yt field name for density is "density"
    #field_name = ('boxlib', 'mag_vel')
    #field_name = ('boxlib', 'y_velocity')
    #field_name = ('boxlib', 'mag_vort')
    field_name = ('boxlib', 'temp')
    #field_name = ('boxlib', 'density')
    #field_name = ('boxlib', 'HeatRelease')
    #field_name = ('boxlib', 'Y(ZMIX)')
    #field_name = ('boxlib', 'MANI_Y-CH4')
    
    
    # Check if this field actually exists in the plotfile
    if field_name not in ds.field_list:
        print(f"  [ERROR] Field '{field_name}' not found in {plt_path}. Skipping.")
        print(f"  Available fields are: {ds.field_list}")
        return

    # Extract data as a fixed-resolution buffer (FRB)
    slc = yt.SlicePlot(ds, "z", field_name)
    frb = slc.frb
    
    # Get coordinates (shift to start at 0,0)
    x = frb["x"].v - ds.domain_left_edge[0].v  # Shift to origin
    y = frb["y"].v - ds.domain_left_edge[1].v  # Shift to origin
    
    # Get the actual data
    data = frb[field_name].v
    
    # Create figure and axes
    fig, ax = plt.subplots(figsize=(8, 8))  # Square figure
    
    img = ax.imshow(
        data,
       # extent=[x.min(), x.max(), y.min(), y.max()],
        extent = get_safe_extent(x, y),
        origin="lower",  # Origin at bottom-left
        cmap="RdBu_r",
        aspect="equal",  # Enforce 1:1 aspect ratio
    )
    
    plt.colorbar(img, ax=ax, label="")
    plt.hlines(y=10.5*0.0762+0.001+1.0*0.0762, xmin=0, xmax=1.6002, colors='cyan', linestyles='--', label='Extraction Line') #line for Fig 4 validation (Mohit 2025)
    plt.hlines(y=10.5*0.0762+0.001+6.0*0.0762, xmin=0, xmax=1.6002, colors='cyan', linestyles='--', label='Extraction Line') #line for Fig 4 validation (Mohit 2025)
    plt.hlines(y=10.5*0.0762+0.001+11.0*0.0762, xmin=0, xmax=1.6002, colors='cyan', linestyles='--', label='Extraction Line') #line for Fig 4 validation (Mohit 2025)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    # Get current time from the dataset
    sim_time = ds.current_time.in_units('s')
    plot_name = os.path.basename(plt_path)
    
    # Set the title with both name and time
    # .4e format is usually best for simulation time as it uses scientific notation
    ax.set_title(f"File: {plot_name}  Time: {sim_time:.4e}")
    
    # Save
    save_path = f"{output_prefix}.png"
    plt.savefig(save_path, bbox_inches="tight", dpi=300)
    
    # CRITICAL: Close the figure to free memory
    #plt.close(fig)
    plt.open(fig)

# --- Main execution block ---
if __name__ == "__main__":
    
    # ***** THIS IS THE FIX *****
    import argparse
    # ***************************
    
    parser = argparse.ArgumentParser(description='Plot density for all timesteps in a directory.')
    parser.add_argument("base_dir", help="Path to the directory containing pltXXXXX folders")
    parser.add_argument("-o", "--output", default=filename, help="Output filename prefix (e.g., 'Rho')")
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
            plot_rho(plt_path, output_prefix_for_func)
        except Exception as e:
            print(f"  [ERROR] Failed to process {plot_name}: {e}")
    
    print("\n--- Processing complete ---")