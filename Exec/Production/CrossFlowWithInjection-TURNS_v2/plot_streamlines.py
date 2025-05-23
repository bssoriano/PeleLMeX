##!/usr/bin/env python3
#import yt
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import ticker
#
#def plot_velocity_streamlines(plt_path, output_file="velocity_streamlines.png"):
#    # Load plotfile
#    ds = yt.load(plt_path)
#    
#    # Extract data on the finest level
#    ad = ds.all_data()
#    x = ad["boxlib", "x"].to("m").v
#    y = ad["boxlib", "y"].to("m").v
#    velx = ad["boxlib", "x_velocity"].v
#    vely = ad["boxlib", "y_velocity"].v
#    vel_mag = np.sqrt(velx**2 + vely**2)
#    
#    # Create grid for streamlines
#    nx, ny = 100, 100  # Resolution of streamline grid
#    xi = np.linspace(x.min(), x.max(), nx)
#    yi = np.linspace(y.min(), y.max(), ny)
#    X, Y = np.meshgrid(xi, yi)
#    
#    # Interpolate velocities to regular grid
#    from scipy.interpolate import griddata
#    Vx = griddata((x, y), velx, (X, Y), method='linear')
#    Vy = griddata((x, y), vely, (X, Y), method='linear')
#    Vmag = griddata((x, y), vel_mag, (X, Y), method='linear')
#    
#    # Create plot
#    plt.figure(figsize=(10, 8))
#    
#    # Background color shows velocity magnitude
#    plt.pcolormesh(X, Y, Vmag, cmap="viridis", shading='auto')
#    cbar = plt.colorbar(label="Velocity Magnitude (m/s)")
#    
#    # Draw streamlines
#    plt.streamplot(X, Y, Vx, Vy, 
#                  color='white', 
#                  linewidth=0.5, 
#                  arrowsize=1.0,
#                  density=2.0)
#    
#    # Formatting
#    plt.title("Velocity Streamlines")
#    plt.xlabel("x (m)")
#    plt.ylabel("y (m)")
#    plt.gca().set_aspect('equal')
#    
#    # Use scientific notation for large/small scales
#    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
#    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
#    
#    plt.savefig(output_file, dpi=300, bbox_inches='tight')
#    print(f"Saved streamline plot to {output_file}")
#
#if __name__ == "__main__":
#    import argparse
#    parser = argparse.ArgumentParser()
#    parser.add_argument("plt_dir", help="Path to pltXXXXX directory")
#    parser.add_argument("-o", "--output", default="velocity_streamlines.png")
#    args = parser.parse_args()
#    plot_velocity_streamlines(args.plt_dir, args.output)
#
#
#
#

#!/usr/bin/env python3
import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib import ticker

def find_last_plt(directory="."):
    """Find the most recent plt folder in the given directory"""
    plt_dirs = sorted(glob.glob(os.path.join(directory, "plt[0-9][0-9][0-9][0-9][0-9]")))
    if not plt_dirs:
        raise FileNotFoundError("No pltXXXXX directories found")
    return plt_dirs[-1]

def plot_velocity_streamlines(plt_path, output_file="output/velocity_streamlines.png"):
    # Load plotfile
    ds = yt.load(plt_path)
    
    # Extract data on the finest level
    ad = ds.all_data()
    x = ad["boxlib", "x"].v
    y = ad["boxlib", "y"].v
    velx = ad["boxlib", "x_velocity"].v
    vely = ad["boxlib", "y_velocity"].v
    vel_mag = np.sqrt(velx**2 + vely**2)
    
    # Create grid for streamlines
    nx, ny = 100, 100  # Resolution
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    X, Y = np.meshgrid(xi, yi)
    
    # Interpolate velocities
    from scipy.interpolate import griddata
    Vx = griddata((x, y), velx, (X, Y), method='linear')
    Vy = griddata((x, y), vely, (X, Y), method='linear')
    Vmag = griddata((x, y), vel_mag, (X, Y), method='linear')
    
    # Create plot
    plt.figure(figsize=(12, 10))
    
    # Background velocity magnitude
    plt.pcolormesh(X, Y, Vmag, cmap="RdBu_r", shading='auto')
    cbar = plt.colorbar(label="Velocity Magnitude (m/s)")
    
    # Streamlines
    plt.streamplot(X, Y, Vx, Vy, 
                  color='black',
                  linewidth=1.0,
                  arrowsize=1.2,
                  density=2.0,
                  arrowstyle='->')
    
    # Formatting
    plt.title(f"Velocity Streamlines (Time = {ds.current_time.to('s'):.2f})")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.gca().set_aspect('equal')
    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved streamline plot of {os.path.basename(plt_path)} to {output_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Plot streamlines from last plt folder')
    parser.add_argument('-d', '--directory', default='.', help='Directory to search for plt folders')
    parser.add_argument('-o', '--output', default='velocity_streamlines.png', help='Output filename')
    args = parser.parse_args()
    
    try:
        last_plt = find_last_plt(args.directory)
        plot_velocity_streamlines(last_plt, args.output)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)