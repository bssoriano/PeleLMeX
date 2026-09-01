import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================================
# 1. Configuration & Parameters
# ==========================================
plotfile = "plt19780"
iso_variable = "temp"      
isovalue = 800.0            

# Burner Cylinder Parameters
cyl_radius = 0.0481
cyl_height = 0.762
cyl_center_x = 0.0    
cyl_center_y = 0.762  
cyl_center_z = 0.0    

# Domain Boundaries
x_lo, x_hi = 0.0, 1.6002
y_lo, y_hi = 0.0, 4.8006
z_lo, z_hi = -0.8001, 0.8001

# ==========================================
# 2. Extract the 3D "Glass Shell" (T = 800 K)
# ==========================================
print(f"Loading {plotfile}...")
ds = yt.load(plotfile)
ds.force_periodicity()

ad = ds.all_data()
print(f"Extracting 3D surface for {iso_variable} = {isovalue} K...")
surface = ds.surface(ad, iso_variable, isovalue)

# Map vertices: Plot X = Sim Z | Plot Y = Sim X | Plot Z = Sim Y
verts = surface.vertices.d
verts_swapped = np.zeros_like(verts)
verts_swapped[0] = verts[2]  # Plot X (Width)  <- Sim Z
verts_swapped[1] = verts[0]  # Plot Y (Depth)  <- Sim X
verts_swapped[2] = verts[1]  # Plot Z (Height) <- Sim Y
triangles = verts_swapped.T.reshape(-1, 3, 3)

# ==========================================
# 3. Extract the 2D Temperature Slice
# ==========================================
print("Extracting 2D Temperature slice at Z = 0.0...")
slc = ds.slice('z', 0.0)

# Create a grid for the slice (higher res for Y axis)
N_x, N_y = 200, 600
# Pass as a tuple of tuples with explicit units to prevent yt from getting confused
width = ( (x_hi - x_lo, 'code_length'), (y_hi - y_lo, 'code_length') )
center = ((x_lo + x_hi) / 2.0, (y_lo + y_hi) / 2.0, 0.0)

frb = slc.to_frb(width, (N_x, N_y), center=center)

# Flip the array vertically so it aligns with standard Matplotlib grid orientation
temp_slice = np.flipud(np.array(frb['temp']))

# Create matching meshgrid for the contour projection
sim_x_1d = np.linspace(x_lo, x_hi, N_x)
sim_y_1d = np.linspace(y_lo, y_hi, N_y)
X_grid, Y_grid = np.meshgrid(sim_x_1d, sim_y_1d)

# ==========================================
# 4. Setup Matplotlib 3D Hybrid Plot
# ==========================================
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# A. Plot the 2D Temperature Slice INSIDE the 3D Domain
# zdir='x' tells Matplotlib to project this onto the Plot X-axis (which is Sim Z = 0.0)
# Therefore, the X_grid (Sim X) maps to Plot Y, and Y_grid (Sim Y) maps to Plot Z
slice_plot = ax.contourf(X_grid, Y_grid, temp_slice, zdir='x', offset=0.0, 
                         levels=60, cmap='turbo', alpha=0.9, antialiased=True)

# Add a colorbar for the internal temperature
cbar = fig.colorbar(slice_plot, ax=ax, shrink=0.5, pad=0.1)
cbar.set_label('Internal Temperature [K]', fontsize=12, weight='bold')

# B. Plot the Transparent 3D "Glass Shell"
# We use a white/light grey color with very low alpha so we can see the slice inside
flame_mesh = Poly3DCollection(triangles, alpha=0.15, facecolors='whitesmoke', linewidths=0)
flame_mesh.set_edgecolor('none')
ax.add_collection3d(flame_mesh)

# ==========================================
# 5. Generate the Black Cylinder (EB Geometry)
# ==========================================
x_cyl_sim = np.linspace(cyl_center_x - cyl_height/2.0, cyl_center_x + cyl_height/2.0, 50)
theta = np.linspace(0, 2 * np.pi, 50)
theta_grid, x_grid_sim = np.meshgrid(theta, x_cyl_sim)

y_grid_sim = cyl_radius * np.cos(theta_grid) + cyl_center_y
z_grid_sim = cyl_radius * np.sin(theta_grid) + cyl_center_z

plot_x_cyl = z_grid_sim  
plot_y_cyl = x_grid_sim  
plot_z_cyl = y_grid_sim  
ax.plot_surface(plot_x_cyl, plot_y_cyl, plot_z_cyl, color='black', alpha=1.0, shade=True)

# ==========================================
# 6. Format Axes and Display
# ==========================================
px_lo, px_hi = z_lo, z_hi  
py_lo, py_hi = x_lo, x_hi  
pz_lo, pz_hi = y_lo, y_hi  

ax.set_xlim(px_lo, px_hi)
ax.set_ylim(py_lo, py_hi)
ax.set_zlim(pz_lo, pz_hi)

# Clean up backgrounds
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.grid(False)

ax.set_xlabel('Z [m] (Width)')
ax.set_ylabel('X [m] (Depth)')
ax.set_zlabel('Y [m] (Height)')
ax.set_title(f"3D Flame Surface (T={isovalue}K) with Internal Temperature Field")

# The perfect side-profile isometric view
ax.view_init(elev=15., azim=75)

try:
    ax.set_box_aspect((px_hi - px_lo, py_hi - py_lo, pz_hi - pz_lo))
except AttributeError:
    pass 

print("Rendering hybrid plot...")
plt.tight_layout()
plt.savefig("hybrid_flame_glass_shell.png", dpi=300, bbox_inches='tight')
plt.show()