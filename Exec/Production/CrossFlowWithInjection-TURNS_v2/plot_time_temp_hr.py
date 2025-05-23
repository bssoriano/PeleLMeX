#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from matplotlib import ticker, patches

def create_output_folder(script_dir):
    """Create output folder if it doesn't exist"""
    output_folder = os.path.join(script_dir, "output")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    return output_folder

def plot_temp_and_heat_release(plt_path, output_folder):
    """Plot Temperature and HeatRelease side-by-side with AMReX mesh overlay"""
    try:
        ds = yt.load(plt_path)
        
        # Create a 1x2 subplot figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # --- Plot Temperature with AMReX Mesh (Left Panel) ---
        slc_temp = yt.SlicePlot(ds, "z", ("boxlib", "temp"))
        frb_temp = slc_temp.frb
        
        x = frb_temp["x"].v - ds.domain_left_edge[0].v
        y = frb_temp["y"].v - ds.domain_left_edge[1].v
        temp = frb_temp[("boxlib", "temp")].v
        
        im1 = ax1.imshow(
            temp,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
        )

        ### Overlay AMReX mesh levels
        for level in range(ds.max_level + 1):
            grids = ds.index.select_grids(level)
            for grid in grids:
                left_edge = grid.LeftEdge[:2].v - ds.domain_left_edge[:2].v
                right_edge = grid.RightEdge[:2].v - ds.domain_left_edge[:2].v
                width = right_edge - left_edge
                rect = patches.Rectangle(
                    left_edge,
                    width[0],
                    width[1],
                    linewidth=0.5,
                    edgecolor='black',
                    facecolor='none',
                    alpha=0.5 - 0.4*level/ds.max_level  # Fainter for finer levels
                )
                ax1.add_patch(rect)

        fig.colorbar(im1, ax=ax1, label="Temperature (K)")
        ax1.set_xlabel("x (m)")
        ax1.set_ylabel("y (m)")
        ax1.set_title("Temperature Field with AMReX Mesh")
        
        # --- Plot HeatRelease (Right Panel) ---
        slc_hr = yt.SlicePlot(ds, "z", ("boxlib", "HeatRelease"))
        frb_hr = slc_hr.frb
        
        hr = frb_hr[("boxlib", "HeatRelease")].v
        
        im2 = ax2.imshow(
            hr,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
        )
        fig.colorbar(im2, ax=ax2, label="Heat Release (W/m³)")
        ax2.set_xlabel("x (m)")
        ax2.set_ylabel("y (m)")
        ax2.set_title("Heat Release Field")
        
        # Save with unique name
        plt_num = os.path.basename(plt_path)
        output_path = os.path.join(output_folder, f"{plt_num}_temp_hr_mesh.png")
        plt.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close()  # Prevent memory leaks
        print(f"Saved plot to {output_path}")

    except Exception as e:
        print(f"Error processing {plt_path}: {str(e)}")

def process_all_plt_files(script_dir):
    """Find and process all plt* folders"""
    output_folder = create_output_folder(script_dir)
    plt_files = sorted(glob.glob(os.path.join(script_dir, "plt[0-9][0-9][0-9][0-9][0-9]")))
    
    if not plt_files:
        print(f"No plt* folders found in: {script_dir}")
        return
    
    print(f"Found {len(plt_files)} plt* folders. Processing...")
    for plt_file in plt_files:
        plot_temp_and_heat_release(plt_file, output_folder)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script's directory
    print(f"Looking for plt* folders in: {script_dir}")
    process_all_plt_files(script_dir)
    print("\nDone. All plots saved to: 'output/'")