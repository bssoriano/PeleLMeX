import yt

# 1. Load the AMReX/PeleLMeX plotfile
plotfile = 'plt00068'
ds = yt.load(plotfile)

# 2. Choose a field to plot as the background. 
# 'density' or 'temp' are almost always present in PeleLMeX.
field = 'temp'

# 3. Create a Slice Plot. 
# 'z' means we are slicing normal to the Z-axis (viewing the X-Y plane). 
# For a 2D simulation, 'z' is exactly what you need.
p = yt.SlicePlot(ds, normal='z', fields=field)

# 4. Annotate the individual computational cells
# This draws the exact mesh lines of every single cell.
# Adjust alpha and line_width as needed so it doesn't obscure the data.
p.annotate_cell_edges(line_width=0.002, alpha=0.5, color='white')

# 5. Optional: Annotate the AMR block boundaries (Grids)
# This draws thicker boxes around the distinct AMR refinement patches,
# which is very useful to see how PeleLMeX clustered the refinement.
p.annotate_grids(alpha=0.8, min_level=0, edgecolors='black', linewidth=1.5)

# 6. Save the resulting image
output_filename = 'amrex_mesh_visualization.png'
p.save(output_filename)

print(f"Plot saved to {output_filename}")