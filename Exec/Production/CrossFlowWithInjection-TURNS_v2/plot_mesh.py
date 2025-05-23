##!/usr/bin/env python3
#import yt
#import numpy as np
#import matplotlib.pyplot as plt
#import glob
#import os
#from matplotlib import ticker, colors
#
#def find_last_plt(directory="."):
#    """Find the most recent plt folder in the given directory"""
#    plt_dirs = sorted(glob.glob(os.path.join(directory, "plt[0-9][0-9][0-9][0-9][0-9]")))
#    if not plt_dirs:
#        raise FileNotFoundError("No pltXXXXX directories found")
#    return plt_dirs[-1]
#
#def plot_velocity_with_grid_points(plt_path, output_prefix="velocity_grid_points"):
#    # Load plotfile
#    ds = yt.load(plt_path)
#    
#    # Extract data on the finest level for velocity
#    ad = ds.all_data()
#    x = ad["boxlib", "x"].to("m").v
#    y = ad["boxlib", "y"].to("m").v
#    velx = ad["boxlib", "x_velocity"].v
#    vely = ad["boxlib", "y_velocity"].v
#
#    # Create grid for velocity interpolation
#    nx, ny = 100, 100
#    xi = np.linspace(x.min(), x.max(), nx)
#    yi = np.linspace(y.min(), y.max(), ny)
#    X, Y = np.meshgrid(xi, yi)
#    
#    # Interpolate velocities
#    from scipy.interpolate import griddata
#    Vx = griddata((x, y), velx, (X, Y), method='linear')
#    Vy = griddata((x, y), vely, (X, Y), method='linear')
#
#    # --- Plot X-Velocity with Grid Points ---
#    fig, ax = plt.subplots(figsize=(12, 10))
#    im = ax.pcolormesh(X, Y, Vx, cmap="coolwarm", shading='auto')
#    cbar = fig.colorbar(im, ax=ax, label="X-Velocity (m/s)")
#    
#    # Collect grid points and levels (FIXED: Handle variable grid sizes)
#    grid_points = []
#    grid_levels = []
#    for level in range(ds.max_level + 1):
#        for grid in ds.index.select_grids(level):
#            # Get cell-centered coordinates (flatten to 1D)
#            x_centers = grid["boxlib", "x"].to("m").v.flatten()
#            y_centers = grid["boxlib", "y"].to("m").v.flatten()
#            grid_points.append(np.column_stack((x_centers, y_centers)))
#            grid_levels.append(np.full(len(x_centers), level))
#    
#    # Check if any grid points were collected (FIXED: Use len() instead of truthiness)
#    if len(grid_points) > 0:
#        grid_points = np.concatenate(grid_points)
#        grid_levels = np.concatenate(grid_levels)
#        
#        # Overlay grid points (colored by level)
#        sc = ax.scatter(
#            grid_points[:, 0],
#            grid_points[:, 1],
#            c=grid_levels,
#            cmap="viridis",
#            s=5,
#            alpha=0.6,
#            label="AMReX Grid Points"
#        )
#        cbar_levels = fig.colorbar(sc, ax=ax, label="Refinement Level")
#        ax.legend()
#    
#    # Formatting
#    ax.set_title(f"X-Velocity with AMReX Grid Points (Time = {ds.current_time.to('s'):.2f} s)")
#    ax.set_xlabel("x (m)")
#    ax.set_ylabel("y (m)")
#    ax.set_aspect("equal")
#    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
#    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.3f}"))
#    
#    plt.savefig(f"{output_prefix}_x.png", dpi=300, bbox_inches="tight")
#    plt.close()
#
#    # --- Plot Y-Velocity with Grid Points ---
#    fig, ax = plt.subplots(figsize=(12, 10))
#    im = ax.pcolormesh(X, Y, Vy, cmap="coolwarm", shading='auto')
#    cbar = fig.colorbar(im, ax=ax, label="Y-Velocity (m/s)")
#    
#    if len(grid_points) > 0:  # Reuse points from previous plot
#        sc = ax.scatter(
#            grid_points[:, 0],
#            grid_points[:, 1],
#            c=grid_levels,
#            cmap="viridis",
#            s=5,
#            alpha=0.6,
#            label="AMReX Grid Points"
#        )
#        cbar_levels = fig.colorbar(sc, ax=ax, label="Refinement Level")
#        ax.legend()
#    
#    ax.set_title(f"Y-Velocity with AMReX Grid Points (Time = {ds.current_time.to('s'):.2f} s)")
#    ax.set_xlabel("x (m)")
#    ax.set_ylabel("y (m)")
#    ax.set_aspect("equal")
#    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
#    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.3f}"))
#    
#    plt.savefig(f"{output_prefix}_y.png", dpi=300, bbox_inches="tight")
#    plt.close()
#
#    print(f"Saved X-velocity + grid points to {output_prefix}_x.png")
#    print(f"Saved Y-velocity + grid points to {output_prefix}_y.png")
#
#if __name__ == "__main__":
#    import argparse
#    parser = argparse.ArgumentParser(description='Plot velocities with AMReX grid points')
#    parser.add_argument('-d', '--directory', default='.', help='Directory to search for plt folders')
#    parser.add_argument('-o', '--output', default='velocity_grid_points', help='Output filename prefix')
#    args = parser.parse_args()
#    
#    try:
#        last_plt = find_last_plt(args.directory)
#        plot_velocity_with_grid_points(last_plt, args.output)
#    except Exception as e:
#        print(f"Error: {e}")
#        exit(1)


########################## PLOT AMREX MESH SECTORS: ##############
#!/usr/bin/env python3
import yt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib import ticker, patches

def find_last_plt(directory="."):
    """Find the most recent plt folder in the given directory"""
    plt_dirs = sorted(glob.glob(os.path.join(directory, "plt[0-9][0-9][0-9][0-9][0-9]")))
    if not plt_dirs:
        raise FileNotFoundError("No pltXXXXX directories found")
    return plt_dirs[-1]

def plot_velocity_with_mesh(plt_path, output_prefix="velocity_with_mesh"):
    # Load plotfile
    ds = yt.load(plt_path)
    
    # Extract data on the finest level
    ad = ds.all_data()
    x = ad["boxlib", "x"].to("m").v  # Ensure units are meters
    y = ad["boxlib", "y"].to("m").v
    velx = ad["boxlib", "x_velocity"].v
    vely = ad["boxlib", "y_velocity"].v

    # Create grid for interpolation
    nx, ny = 100, 100  # Resolution
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    X, Y = np.meshgrid(xi, yi)
    
    # Interpolate velocities
    from scipy.interpolate import griddata
    Vx = griddata((x, y), velx, (X, Y), method='linear')
    Vy = griddata((x, y), vely, (X, Y), method='linear')

    # Get AMReX grid patches (for mesh overlay)
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Plot X-Velocity heatmap
    im = ax.pcolormesh(X, Y, Vx, cmap="coolwarm", shading='auto')
    cbar = fig.colorbar(im, ax=ax, label="X-Velocity (m/s)")
    
    # Overlay AMReX mesh
    for grid in ds.index.grids:
        left_edge = grid.LeftEdge.to("m").v
        right_edge = grid.RightEdge.to("m").v
        width = right_edge[0] - left_edge[0]
        height = right_edge[1] - left_edge[1]
        rect = patches.Rectangle(
            (left_edge[0], left_edge[1]),
            width,
            height,
            linewidth=0.5,
            edgecolor="black",
            facecolor="none",
            alpha=0.3,
        )
        ax.add_patch(rect)
    
    # Formatting
    ax.set_title(f"X-Velocity with AMReX Mesh (Time = {ds.current_time.to('s'):.2f} s)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.3f}"))
    
    plt.savefig(f"{output_prefix}_x.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Repeat for Y-Velocity
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.pcolormesh(X, Y, Vy, cmap="coolwarm", shading='auto')
    cbar = fig.colorbar(im, ax=ax, label="Y-Velocity (m/s)")
    
    # Overlay AMReX mesh
    for grid in ds.index.grids:
        left_edge = grid.LeftEdge.to("m").v
        right_edge = grid.RightEdge.to("m").v
        width = right_edge[0] - left_edge[0]
        height = right_edge[1] - left_edge[1]
        rect = patches.Rectangle(
            (left_edge[0], left_edge[1]),
            width,
            height,
            linewidth=0.5,
            edgecolor="black",
            facecolor="none",
            alpha=0.3,
        )
        ax.add_patch(rect)
    
    ax.set_title(f"Y-Velocity with AMReX Mesh (Time = {ds.current_time.to('s'):.2f} s)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.3f}"))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{y:.3f}"))
    
    plt.savefig(f"{output_prefix}_y.png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved X-velocity + mesh plot to {output_prefix}_x.png")
    print(f"Saved Y-velocity + mesh plot to {output_prefix}_y.png")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Plot velocity components with AMReX mesh overlay')
    parser.add_argument('-d', '--directory', default='.', help='Directory to search for plt folders')
    parser.add_argument('-o', '--output', default='velocity_with_mesh', help='Output filename prefix')
    args = parser.parse_args()
    
    try:
        last_plt = find_last_plt(args.directory)
        plot_velocity_with_mesh(last_plt, args.output)
    except Exception as e:
        print(f"Error: {e}")
        exit(1)