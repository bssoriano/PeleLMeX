import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================================
# 1. Configuration & Parameters
# ==========================================
plotfile = "plt19780"
flame_variable = "temp"      # Try "HeatRelease", "temp", or "MANI_Y-OH"
isovalue = 1500.0            # Set to your desired threshold

# Burner Cylinder Parameters (Adjust to match your domain)
cyl_radius = 0.02            # Radius of the burner in meters
cyl_height = 0.05            # Height of the visible black base in meters
cyl_center_x = 0.0           # X coordinate of the cylinder center
cyl_center_y = 0.0           # Y coordinate of the cylinder center

# ==========================================
# 2. Load Data and Extract Flame Surface
# ==========================================
print(f"Loading {plotfile}...")
ds = yt.load(plotfile)

# FIX: Force periodicity to avoid boundary ghost cell read errors
ds.force_periodicity()

ad = ds.all_data()
print(f"Extracting {flame_variable} surface at {isovalue}...")
surface = ds.surface(ad, flame_variable, isovalue)

# ---------------------------------------------------------
# THIS IS WHAT WENT MISSING: Extracting the triangles
# ---------------------------------------------------------
verts = surface.vertices.d
triangles = verts.T.reshape(-1, 3, 3)
print(f"Extracted {len(triangles)} triangles.")

# ==========================================
# 3. Setup Matplotlib 3D Plot
# ==========================================
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Add the Flame Surface
flame_mesh = Poly3DCollection(triangles, alpha=0.6, facecolors='darkorange', linewidths=0)
ax.add_collection3d(flame_mesh)

# ==========================================
# 4. Generate the Black Cylinder (Base)
# ==========================================
z_min = ds.domain_left_edge[2].d
z_max_cyl = z_min + cyl_height

z_cyl = np.linspace(z_min, z_max_cyl, 50)
theta = np.linspace(0, 2 * np.pi, 50)
theta_grid, z_grid = np.meshgrid(theta, z_cyl)

x_grid = cyl_radius * np.cos(theta_grid) + cyl_center_x
y_grid = cyl_radius * np.sin(theta_grid) + cyl_center_y

# Plot the cylinder surface
ax.plot_surface(x_grid, y_grid, z_grid, color='black', alpha=1.0, shade=True)

# ==========================================
# 5. Format Axes and Display
# ==========================================
ax.set_xlim(ds.domain_left_edge[0].d, ds.domain_right_edge[0].d)
ax.set_ylim(ds.domain_left_edge[1].d, ds.domain_right_edge[1].d)
ax.set_zlim(ds.domain_left_edge[2].d, ds.domain_right_edge[2].d)

ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
ax.set_title(f"Flame Isocontour ({flame_variable} = {isovalue})")

ax.view_init(elev=20., azim=45)

print("Rendering plot...")
plt.tight_layout()
plt.show()