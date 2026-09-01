#!/usr/bin/env python3
import yt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import numpy as np
import argparse
import os

def plot_detailed_mesh(plt_path):
    if not os.path.exists(plt_path):
        print(f"Error: {plt_path} not found")
        return

    print(f"Loading data from {plt_path}...")
    ds = yt.load(plt_path)

    # Setup the figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Set plot limits based on the whole domain
    domain_left = ds.domain_left_edge.v
    domain_right = ds.domain_right_edge.v
    
    # Adjust indices if your simulation is not XY (0,1). 
    # For standard AMReX 2D or 3D slice, this is usually 0 (x) and 1 (y).
    ix, iy = 0, 1
    
    ax.set_xlim(domain_left[ix], domain_right[ix])
    ax.set_ylim(domain_left[iy], domain_right[iy])
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title(f"Mesh: {os.path.basename(plt_path)}\nRed = Grid Boundaries | Black = Cells")

    # Storage for the black cell lines (using a list is faster than adding to plot one by one)
    cell_segments = []

    print(f"Processing {ds.index.num_grids} grids...")
    
    # Loop over every AMR grid box
    for g in ds.index.grids:
        # 1. GET COORDINATES
        le = g.LeftEdge.v  # Left edge [x, y, z]
        re = g.RightEdge.v # Right edge [x, y, z]
        dims = g.ActiveDimensions # Number of cells [nx, ny, nz]
        
        # Grid Box Dimensions
        g_width = re[ix] - le[ix]
        g_height = re[iy] - le[iy]

        # 2. DRAW THE RED BOX (The AMR Patch)
        # zorder=10 ensures the red box sits ON TOP of the black lines
        rect = patches.Rectangle((le[ix], le[iy]), g_width, g_height, 
                                 linewidth=2.0, edgecolor='red', facecolor='none', zorder=10)
        ax.add_patch(rect)

        # 3. CALCULATE BLACK CELL LINES
        # We generate the start/end points for the internal grid lines
        
        # Vertical lines (constant X, varying Y)
        # We skip the first and last line because the Red Box already covers the border
        x_steps = np.linspace(le[ix], re[ix], dims[ix] + 1)
        for x in x_steps[1:-1]: 
            cell_segments.append([(x, le[iy]), (x, re[iy])])

        # Horizontal lines (constant Y, varying X)
        y_steps = np.linspace(le[iy], re[iy], dims[iy] + 1)
        for y in y_steps[1:-1]:
            cell_segments.append([(le[ix], y), (re[ix], y)])

    print(f"Rendering {len(cell_segments)} cell lines... (this might pause for a second)")
    
    # 4. ADD ALL BLACK LINES AT ONCE (High Performance)
    lc = LineCollection(cell_segments, colors='black', linewidths=0.5, alpha=0.5, zorder=1)
    ax.add_collection(lc)

    print("Opening Window. Zoom in to see the cells!")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("plt_file", help="Path to plt folder")
    args = parser.parse_args()
    
    plot_detailed_mesh(args.plt_file)