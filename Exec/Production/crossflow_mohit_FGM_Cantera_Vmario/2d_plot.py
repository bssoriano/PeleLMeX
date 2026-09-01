import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ==========================================
# 1. Configuration & Parameters
# ==========================================
plotfile = "plt19780"
background_var = "MANI_Y-CH4"       # The filled color map variable
contour_var = "MANI_Y-OH"     # The overlaid line contour variable
contour_val = 0.0015          # Isovalue for the line contour

# Burner Cylinder Parameters (EB Setup)
cyl_radius = 0.0481
cyl_height = 0.762
cyl_center_x = 0.0    # Sim X
cyl_center_y = 0.762  # Sim Y

# Domain Boundaries
x_lo, x_hi = 0.0, 1.6002
y_lo, y_hi = 0.0, 4.8006
z_lo, z_hi = -0.8001, 0.8001

# ==========================================
# 2. Load Data and Extract 2D Slice
# ==========================================
print(f"Loading {plotfile}...")
ds = yt.load(plotfile)
ds.force_periodicity()

# Cut a slice through the center of the Z-axis (Z=0.0)
# Since the flare is along X and wind along Y, this gives us the X-Y center plane.
print("Slicing data at Z=0.0...")
slc = ds.slice('z', 0.0)

# Create a Fixed Resolution Buffer (FRB) to map the AMR data to a flat 2D grid for Matplotlib
# We use a high resolution (e.g., 1000x1000) for a crisp publication plot
res = (1000, 1000)
# Pass as a tuple of tuples with explicit units to prevent yt from getting confused
width = ((x_hi - x_lo, 'code_length'), (y_hi - y_lo, 'code_length'))
frb = slc.to_frb(width, res, center=( (x_hi+x_lo)/2.0, (y_hi+y_lo)/2.0, 0.0 ))

# Extract the 2D arrays
bg_field = np.array(frb[background_var])
line_field = np.array(frb[contour_var])

# ==========================================
# 3. Setup Matplotlib 2D Plot
# ==========================================
fig, ax = plt.subplots(figsize=(5, 10)) # Taller figure for your long Y-domain

# Define the physical extent of the bounding box so the axes have real meter units
extent = [x_lo, x_hi, y_lo, y_hi]

# Plot the background filled contour (Temperature)
# 'turbo' is the modern, colorblind-friendly replacement for 'jet'
im = ax.imshow(bg_field, origin='lower', extent=extent, cmap='RdBu_r', aspect='equal')
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Temperature [K]', fontsize=12)

# Plot the overlaid flame front contour (OH)
# We make it white so it pops against the colorful background
#ax.contour(line_field, levels=[contour_val], origin='lower', extent=extent, colors='white', linewidths=2.0)

# ==========================================
# 4. Draw the Burner Nozzle
# ==========================================
# In a 2D center-cut (Z=0), the cylinder looks like a rectangle
# It spans from (cyl_center_x - radius) to (cyl_center_x + radius) in X
# And it sits at the Y-boundary (y_lo) and goes up to cyl_center_y
nozzle_width = 0.38 # Because cylinder is oriented along X
nozzle_x_start = cyl_center_x #- (cyl_height/2.0)
nozzle_y_start = 0.762-0.0481#y_lo
nozzle_length = 2*0.0481 # Height of the nozzle sticking into the domain

nozzle = patches.Rectangle((nozzle_x_start, nozzle_y_start), nozzle_width, nozzle_length, 
                           linewidth=1, edgecolor='black', facecolor='black', zorder=10)
ax.add_patch(nozzle)

# ==========================================
# 5. Format Axes and Display
# ==========================================
ax.set_xlabel('$X$ [m]', fontsize=12)
ax.set_ylabel('$Y$ [m]', fontsize=12)
ax.tick_params(axis='both', labelsize=12)

# Set the limits slightly zoomed in if you want to focus on the flame
# ax.set_xlim(0, 1.0)
# ax.set_ylim(0, 3.0)

#plt.title(f"Center Plane Slice (Z=0)\n{background_var} with {contour_var} contour", fontsize=16)

print("Saving 2D publication plot...")
plt.tight_layout()
plt.savefig("flare_2d_slice.png", dpi=300, bbox_inches='tight')
plt.show()