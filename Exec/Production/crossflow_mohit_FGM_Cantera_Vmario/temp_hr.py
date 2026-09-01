#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from matplotlib import ticker, patches
from matplotlib.ticker import MultipleLocator
import numpy.ma as ma
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm

plt.rcParams.update({
    "text.usetex": True,   # Enable LaTeX rendering
    "font.family": "serif", # Use serif font (like LaTeX default)
    "font.serif": ["Arial"],  # Optional: specify font
})

def create_output_folder(script_dir):
    """Create output folder if it doesn't exist"""
    output_folder = os.path.join(script_dir, "output")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    return output_folder

def plot_temp_and_heat_release(plt_path, output_folder):
    """Plot Temperature and HeatRelease as separate figures"""
    try:
        plt_num = os.path.basename(plt_path)
        ds = yt.load(plt_path)
        
        # --- Temperature Plot (Standalone Figure) ---
        fig_temp = plt.figure(figsize=(8, 8))
        ax_temp = fig_temp.gca()
        
        slc_temp = yt.SlicePlot(ds, "z", ("boxlib", "temp"))
        frb_temp = slc_temp.frb
        
        x = frb_temp["x"].v - ds.domain_left_edge[0].v
        y = frb_temp["y"].v - ds.domain_left_edge[1].v
        temp = frb_temp[("boxlib", "temp")].v
        
        # Mask values below 298K
        masked_temp = ma.masked_where(temp < 296, temp)
        
        
        im_temp = ax_temp.imshow(
            masked_temp,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
            vmin=296,
            vmax=temp.max()
        )
        
        # Overlay AMReX mesh levels
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
                    linewidth=0.8,
                    edgecolor='black',
                    facecolor='none',
                    alpha=0.5 - 0.4*level/ds.max_level
                )
                ax_temp.add_patch(rect)

        cbar_temp = fig_temp.colorbar(im_temp)
        cbar_temp.set_label("Temperature (K)")
        tick_step = (temp.max() - 298) / 8
        cbar_temp.set_ticks(np.arange(298, temp.max() + tick_step, tick_step))
        ax_temp.set_xlabel("x (m)")
        ax_temp.set_ylabel("y (m)")
        ax_temp.set_xticks(np.arange(0., x.max(), 0.2))
        ax_temp.set_yticks(np.arange(0., y.max(), 0.4))

        # Define rectangle parameters
        rect_x = 0  # X-coordinate of bottom-left corner (adjust as needed)
        rect_y = 0.774-0.0481  # Y-coordinate of bottom-left corner (adjust as needed)
        rect_width = 0.375  # Length of the rectangle
        rect_height = 2*0.0481  # Height of the rectangle
        print(f"X-axis limits: {plt.xlim()}, Y-axis limits: {plt.ylim()}") 
        # Create the rectangle patch
        rect = Rectangle(
            (rect_x, rect_y),  # Bottom-left corner
            rect_width,        # Width
            rect_height,       # Height
            linewidth=1,       # Border thickness
            edgecolor="black",  # Border color
            facecolor="black",  # Fill color
        )
        ax_temp.add_patch(rect)
        
        # Save Temperature plot
        temp_output = os.path.join(output_folder, f"{plt_num}_temperature_mesh.png")
        fig_temp.savefig(temp_output, bbox_inches="tight", dpi=300)
        plt.close(fig_temp)
        print(f"Saved Temperature plot to {temp_output}")

        # --- Heat Release Plot (Standalone Figure) ---
        fig_hr = plt.figure(figsize=(8, 8))
        ax_hr = fig_hr.gca()
        
        slc_hr = yt.SlicePlot(ds, "z", ("boxlib", "HeatRelease"))
        frb_hr = slc_hr.frb
        
        y_ch4 = frb_hr[("boxlib", "MANI_Y-CH4")].v
        y_o2 = frb_hr[("boxlib", "MANI_Y-O2")].v
        epsilon = 1e-10
        phi = 4.0 * y_ch4 / (y_o2 + epsilon)
        hr = frb_hr[("boxlib", "HeatRelease")].v

        # Ensure hr has no zeros/negatives (add small offset if needed)
        hr_positive = np.maximum(hr, 1e-6)  # Replace values ≤0 with 1e-6
        
        # Set vmin/vmax for log scale (must be >0)
        vmin_log, vmax_log = 1.0, hr_positive.max()  # Adjust based on your data range

        vmin, vmax = 0.0, 1e7  # Set to -6e7 and +6e7
        im_hr = ax_hr.imshow(
            hr_positive,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
            norm=LogNorm(vmin=vmin_log, vmax=vmax_log),  # Key: Apply logarithmic scaling
        )
        
        
        cs = ax_hr.contour(
            x, y, phi, 
            levels=[0.6],
            colors='lime',
            linewidths=1,
            linestyles='--'
        )
        #ax_hr.clabel(cs, cs.levels, inline=True, fmt='CH4/O2 = %1.1f', fontsize=10)
        
        fig_hr.colorbar(im_hr, label="Heat Release (W/m³)")
        # Define rectangle parameters
        rect_x = 0  # X-coordinate of bottom-left corner (adjust as needed)
        rect_y = 0.762-0.0481  # Y-coordinate of bottom-left corner (adjust as needed)
        rect_width = 0.38  # Length of the rectangle
        rect_height = 0.0481  # Height of the rectangle
        print(f"X-axis limits: {plt.xlim()}, Y-axis limits: {plt.ylim()}") 
        # Create the rectangle patch
        rect = Rectangle(
            (rect_x, rect_y),  # Bottom-left corner
            rect_width,        # Width
            rect_height,       # Height
            linewidth=1,       # Border thickness
            edgecolor="black",  # Border color
            facecolor="black",  # Fill color
        )
        ax_hr.add_patch(rect)

        ax_hr.set_xlabel("x (m)")
        ax_hr.set_ylabel("y (m)")
        ax_hr.set_xticks(np.arange(0.0, x.max(), 0.4))
        ax_hr.set_yticks(np.arange(0.0, y.max(), 0.4))
        
        # Save Heat Release plot
        hr_output = os.path.join(output_folder, f"{plt_num}_hr_phi_LFL.png")
        fig_hr.savefig(hr_output, bbox_inches="tight", dpi=300)
        plt.close(fig_hr)
        print(f"Saved Heat Release plot to {hr_output}")

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
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"Looking for plt* folders in: {script_dir}")
    process_all_plt_files(script_dir)
    print("\nDone. All plots saved to: 'output/'")