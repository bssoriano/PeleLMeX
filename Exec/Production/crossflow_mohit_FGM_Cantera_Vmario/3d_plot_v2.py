import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ==========================================
# 1. Configuration & Parameters
# ==========================================
plotfile = "plt19780"
flame_variable = "temp"      
isovalue = 1500.0

# Burner Cylinder Parameters
cyl_radius = 0.02            # Radius of the burner in meters
cyl_height = 0.05            # Height of the visible black base in meters
cyl_center_x = 0.0           # Sim X coordinate of the cylinder center
cyl_center_z = 0.0           # Sim Z coordinate of the cylinder center

# Domain Boundaries from your inputs
x_lo, x_hi = 0.0, 1.6002
y_lo, y_hi = 0.0, 4.8006
z_lo, z_hi = -0.8001, 0.8001

# ==========================================
# 2. Load Data and Extract Flame Surface
# ==========================================
print(f"Loading {plotfile}...")
ds = yt.load(plotfile)

ad = ds.all_data()

# Find the absolute maximum OH mass fraction in the domain
min_val, max_val = ad.quantities.extrema("HeatRelease")
print(f"Maximum HeatRelease in the domain: {max_val.d:.12f}")

ds.force_periodicity()

ad = ds.all_data()
print(f"Extracting {flame_variable} surface at {isovalue}...")
surface = ds.surface(ad, flame_variable, isovalue)

# Extract vertices and apply coordinate mapping for Matplotlib
verts = surface.vertices.d

# NEW MAPPING: 
# Plot X (Width)  = Sim Z
# Plot Y (Depth)  = Sim X
# Plot Z (Height) = Sim Y
verts_swapped = np.zeros_like(verts)
verts_swapped[0] = verts[2]  # Plot X <- Sim Z
verts_swapped[1] = verts[0]  # Plot Y <- Sim X
verts_swapped[2] = verts[1]  # Plot Z <- Sim Y

triangles = verts_swapped.T.reshape(-1, 3, 3)
print(f"Extracted {len(triangles)} triangles.")

# ==========================================
# 3. Setup Matplotlib 3D Plot
# ==========================================
fig = plt.figure(figsize=(8,16))
ax = fig.add_subplot(111, projection='3d')

# Add the Flame Surface
flame_mesh = Poly3DCollection(triangles, alpha=0.7, facecolors='darkorange', linewidths=0)
#flame_mesh = Poly3DCollection(triangles, alpha=0.5, facecolors='dodgerblue', linewidths=0)
ax.add_collection3d(flame_mesh)

# ==========================================
# 4. Generate the Black Cylinder (EB Geometry)
# ==========================================
# Parameters from PeleLMex EB Setup
cyl_radius = 0.0481
cyl_height = 0.762
cyl_center_x = 0.0    # Sim X
cyl_center_y = 0.762  # Sim Y
cyl_center_z = 0.0    # Sim Z

# Since direction = 0, the cylinder is aligned along Simulation X.
# AMReX centers the cylinder exactly at cyl_center_x.
x_cyl_sim = np.linspace(cyl_center_x, 
                        cyl_center_x + cyl_height/2.0, 50)

theta = np.linspace(0, 2 * np.pi, 50)
theta_grid, x_grid_sim = np.meshgrid(theta, x_cyl_sim)

# The circular cross-section is in the Sim Y and Sim Z plane
y_grid_sim = cyl_radius * np.cos(theta_grid) + cyl_center_y
z_grid_sim = cyl_radius * np.sin(theta_grid) + cyl_center_z

# Map cylinder coordinates to Plot coordinates
# Plot X (Width)  = Sim Z
# Plot Y (Depth)  = Sim X
# Plot Z (Height) = Sim Y
plot_x_cyl = z_grid_sim  
plot_y_cyl = x_grid_sim  
plot_z_cyl = y_grid_sim  

# Plot the cylinder surface
ax.plot_surface(plot_x_cyl, plot_y_cyl, plot_z_cyl, color='black', alpha=1.0, shade=True)


# ==========================================
# 5. Draw the Red Domain Bounding Box
# ==========================================
# Map the domain bounds to the plot axes
px_lo, px_hi = z_lo, z_hi  # Plot X = Sim Z
py_lo, py_hi = x_lo, x_hi  # Plot Y = Sim X
pz_lo, pz_hi = y_lo, y_hi  # Plot Z = Sim Y

# Bottom rectangle (at min Plot Z)
ax.plot([px_lo, px_hi, px_hi, px_lo, px_lo], 
        [py_lo, py_lo, py_hi, py_hi, py_lo], 
        [pz_lo, pz_lo, pz_lo, pz_lo, pz_lo], color='dimgrey', linewidth=1.5, linestyle='--')

# Top rectangle (at max Plot Z)
ax.plot([px_lo, px_hi, px_hi, px_lo, px_lo], 
        [py_lo, py_lo, py_hi, py_hi, py_lo], 
        [pz_hi, pz_hi, pz_hi, pz_hi, pz_hi], color='dimgrey', linewidth=1.5, linestyle='--')

# 4 Vertical pillars connecting top and bottom
for px in [px_lo, px_hi]:
    for py in [py_lo, py_hi]:
        ax.plot([px, px], [py, py], [pz_lo, pz_hi], color='dimgrey', linewidth=1.5, linestyle='--')

# ==========================================
# 6. Format Axes and Display
# ==========================================
# Set axes limits
ax.set_xlim(px_lo, px_hi)
ax.set_ylim(py_lo, py_hi)  # If you still want Depth inverted, change to: ax.set_ylim(py_hi, py_lo)
ax.set_zlim(pz_lo, pz_hi)

# Keep the physical aspect ratio realistic based on the new mapping
try:
    ax.set_box_aspect((px_hi - px_lo, py_hi - py_lo, pz_hi - pz_lo))
except AttributeError:
    pass 

# Label the axes based on the mapped simulation coordinate names
ax.set_xlabel('Z [m]')
ax.set_ylabel('X [m]')
ax.set_zlabel('Y [m]')
ax.set_title(f"Flame Isocontour (T=1500K)")

ax.set_yticks(np.arange(0.0, 1.61, 0.4))
ax.set_xticks(np.arange(-0.8, 0.8, 0.4))

ax.set_aspect('equal')

# 4. Clean up the axes background panes (makes it look much cleaner)
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.grid(False) # Turn off the heavy grid lines

# Adjust view angle 
ax.view_init(elev=20., azim=35)

print("Rendering plot...")
# 6. Save in High-Resolution for the Report!
print("Saving publication-ready plot...")
plt.tight_layout()
plt.savefig("flare_temp_isocontour.png", dpi=300, bbox_inches='tight', transparent=False)
plt.show()

