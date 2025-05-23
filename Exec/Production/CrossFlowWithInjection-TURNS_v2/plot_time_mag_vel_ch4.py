#!/usr/bin/env python3
import yt
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from matplotlib import ticker, patches
from scipy.interpolate import griddata

def create_output_folder(script_dir):
    """Create output folder if it doesn't exist"""
    output_folder = os.path.join(script_dir, "output")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    return output_folder

def plot_mag_vel_and_heat_release(plt_path, output_folder):
    """Plot Velocity Magnitude with Streamlines and CH4 side-by-side"""
    try:
        ds = yt.load(plt_path)
        
        # Create a 1x2 subplot figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # --- Left Panel: Velocity Magnitude with Streamlines ---
        # Extract all data (not just slice)
        ad = ds.all_data()
        x = ad["boxlib", "x"].v
        y = ad["boxlib", "y"].v
        velx = ad["boxlib", "x_velocity"].v
        vely = ad["boxlib", "y_velocity"].v
        vel_mag = np.sqrt(velx**2 + vely**2)
        
        # Create interpolation grid
        nx, ny = 100, 100
        xi = np.linspace(x.min(), x.max(), nx)
        yi = np.linspace(y.min(), y.max(), ny)
        X, Y = np.meshgrid(xi, yi)
        
        # Interpolate velocities
        Vx = griddata((x, y), velx, (X, Y), method='linear')
        Vy = griddata((x, y), vely, (X, Y), method='linear')
        Vmag = griddata((x, y), vel_mag, (X, Y), method='linear')
        
        # Plot velocity magnitude
        im1 = ax1.pcolormesh(X, Y, Vmag, cmap="RdBu_r", shading='auto')
        fig.colorbar(im1, ax=ax1, label="Velocity Magnitude (m/s)")
        
        # Add streamlines
        ax1.streamplot(
            X, Y, Vx, Vy,
            color='black',
            linewidth=0.8,
            arrowsize=1.2,
            density=2.0,
            arrowstyle='->'
        )
        ax1.set_title(f"Velocity (t = {ds.current_time.to('s'):.2f} s)")
        ax1.set_xlabel("x (m)")
        ax1.set_ylabel("y (m)")
        ax1.set_aspect('equal')
        
        # --- Right Panel: CH4 Mass Fraction ---
        # Get slice data for CH4
        slc = yt.SlicePlot(ds, "z", ("boxlib", "Y(CH4)"))
        frb = slc.frb
        x_ch4 = frb["x"].v - ds.domain_left_edge[0].v
        y_ch4 = frb["y"].v - ds.domain_left_edge[1].v
        ch4 = frb[("boxlib", "Y(CH4)")].v
        
        im2 = ax2.imshow(
            ch4,
            extent=[x_ch4.min(), x_ch4.max(), y_ch4.min(), y_ch4.max()],
            origin="lower",
            cmap="RdBu_r",
            aspect="equal",
        )
        fig.colorbar(im2, ax=ax2, label="CH4 Mass Fraction")
        ax2.set_xlabel("x (m)")
        ax2.set_ylabel("y (m)")
        ax2.set_title("CH4 Distribution")
        
        # Save plot
        plt_num = os.path.basename(plt_path)
        output_path = os.path.join(output_folder, f"{plt_num}_velocity_ch4.png")
        plt.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close()
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
        plot_mag_vel_and_heat_release(plt_file, output_folder)

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"Looking for plt* folders in: {script_dir}")
    process_all_plt_files(script_dir)
    print("\nDone. All plots saved to: 'output/'")