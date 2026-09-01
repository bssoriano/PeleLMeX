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
    try:
        plt_num = os.path.basename(plt_path)
        ds = yt.load(plt_path)
        
        

        # --- Heat Release Plot (Standalone Figure) ---
        fig = plt.figure(figsize=(8, 8))
        ax = fig.gca()
        
        slc = yt.SlicePlot(ds, "z", ("boxlib", "HeatRelease"))
        frb = slc.frb
        
        x = frb["x"].v - ds.domain_left_edge[0].v
        y = frb["y"].v - ds.domain_left_edge[1].v
        y_ch4 = frb[("boxlib", "MANI_Y-CH4")].v
        y_o2 = frb[("boxlib", "MANI_Y-O2")].v
        y_n2 = frb[("boxlib", "MANI_Y-N2")].v
        
        # from turns:
        # phi = (A/F)_stoic / (A/F), A/F is the air fuel ratio

        AF_stoic = 17.11 #for CH4 - air
        AF = (y_o2 + y_n2)/y_ch4       
        phi = AF_stoic/AF


        #epsilon = 1e-10
        #phi = 4.0 * y_ch4 / (y_o2 + epsilon)

        #if y_ch4 > 0.0044: #LFL = 4.4%
        #    return 0.0  # Too rich to burn (above LFL)
        #else:
        #    return y_ch4  # Below or at LFL (flammable)
       
        phi_LFL = 0.6  # Lower Flammability Limit (4.4%)
        # Calculate equivalence ratio (phi) from y_ch4
        # Assumes y_air = 1 - y_ch4 (for simplicity; adjust if needed)
        #stoich_ratio = 0.055  # Stoichiometric CH4/air ratio (~9.5% CH4 in air)
        #phi = (y_ch4 / (1 - y_ch4)) / stoich_ratio
        
        # Mask unsafe regions (phi >= phi_LFL or y_ch4 > LFL)
        mask = (phi < phi_LFL) & (y_ch4 > 0)  # Exclude phi >= phi_LFL and zero/negative values
        y_ch4_safe = np.where(mask, y_ch4, 0.0)  # Set unsafe regions to 0.0 

        #y_ch4_safe = np.where(phi_LFL > phi, 0.0, y_ch4)  # Set >LFL to 0

        vmin_log, vmax_log = 1.e-6, y_ch4_safe.max()  # Adjust based on your data range
        im = ax.imshow(
            y_ch4_safe,
            extent=[x.min(), x.max(), y.min(), y.max()],
            origin="lower",
            cmap="OrRd",
            aspect="equal",
            norm=LogNorm(vmin=vmin_log, vmax=vmax_log)  # Key: Apply logarithmic scaling
        )
        
        
        cs = ax.contour(
            x, y, phi, 
            levels=[phi_LFL],
            colors='black',
            linewidths=0.6,
            linestyles='-'
        )
        
        fig.colorbar(im, label="$Y_{CH_4}$")
        #fmt = '%1.4f'
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
        ax.add_patch(rect)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_xticks(np.arange(0.0, x.max(), 0.2))
        ax.set_yticks(np.arange(0.0, y.max(), 0.4))
        #ax.set_xlim(0,0.4)
        #ax.set_ylim(0,0.6)
        
        # Save Heat Release plot
        hr_output = os.path.join(output_folder, f"{plt_num}_ch4_phi.png")
        fig.savefig(hr_output, bbox_inches="tight", dpi=300)
        plt.close(fig)
        print(f"Saved Heat Release plot to {hr_output}")

    except Exception as e:
        print(f"Error processing {plt_path}: {str(e)}")

def process_all_plt_files(script_dir):
    """Find and process all plt* folders"""
    output_folder = create_output_folder(script_dir)
    plt_files = sorted(glob.glob(os.path.join(script_dir, "plt[0-9][0-9][0-9][0-9][0-9]*")))
    
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