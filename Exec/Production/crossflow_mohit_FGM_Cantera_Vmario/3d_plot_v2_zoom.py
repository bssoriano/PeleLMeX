import numpy as np
import pyvista as pv
import yt

# ==============================================================================
# 1. PARAMETERS & INPUT CONFIGURATION
# ==============================================================================
# Path to your PeleLMeX plotfile directory
PLT_DIR = "plt41967"

# Target field parameter in BoxLib format for PeleLMeX
TEMP_FIELD = ("boxlib", "temp")

# Flame reaction zone threshold (Kelvin)
T_ISO_VALUE = 800.0
T_AMBIENT = 300.0  # Floor for color mapping normalization

# Burner Embedded Boundary (EB2) Parameters from PeleLMeX input:
EB_RADIUS = 0.0481
EB_HEIGHT = 0.762
EB_CENTER = (0.0, 0.762, 0.0)
EB_DIRECTION_AXIS = 0  # X-axis

# ==============================================================================
# 2. BOXLIB DATA EXTRACTION VIA YT
# ==============================================================================
print(f"[INFO] Loading PeleLMeX dataset: {PLT_DIR}")
ds = yt.load(PLT_DIR)

# Extract domain boundaries and grid dimensions at max refinement level
domain_left = ds.domain_left_edge.value  # SI Units (m)
domain_right = ds.domain_right_edge.value
max_level = ds.max_level
ref_dims = ds.domain_dimensions * (2**max_level)

print(
    f"[INFO] Domain Bounds: {domain_left} to {domain_right} | Max AMR Level: {max_level}"
)

# Extract dynamic covering grid at max refinement level
cg = ds.covering_grid(level=max_level, left_edge=domain_left, dims=ref_dims)

# Extract 3D scalar temperature array (Kelvin)
temp_array = cg[TEMP_FIELD].v
temp_max = float(temp_array.max())
print(f"[INFO] Temperature field extracted: Min = {temp_array.min():.1f} K, Max = {temp_max:.1f} K")

# ==============================================================================
# 3. PYVISTA MESH, THRESHOLDING & CLIPPING
# ==============================================================================
print("[INFO] Constructing PyVista ImageData structure...")
grid = pv.ImageData()
grid.dimensions = np.array(temp_array.shape) + 1
grid.origin = domain_left
grid.spacing = (domain_right - domain_left) / temp_array.shape

# Assign cell scalar data and convert to point data for smooth spatial interpolation
grid.cell_data["temp"] = temp_array.flatten(order="F")
grid_points = grid.cell_data_to_point_data()

# 1. Extract 3D internal flame volume where T >= 800 K
print(f"[INFO] Thresholding grid for internal volume T >= {T_ISO_VALUE} K...")
flame_volume = grid_points.threshold(value=T_ISO_VALUE, scalars="temp")

# 2. Clip thresholded volume along the z = 0 centerplane to expose interior temperature
print("[INFO] Applying centerplane cut at z = 0 m...")
clipped_volume = flame_volume.clip(
    normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 0.0)
)

# 3. Generate outer 800 K boundary shell for reference
print(f"[INFO] Generating T = {T_ISO_VALUE} K boundary shell...")
flame_shell = grid_points.contour(
    scalars="temp", isosurfaces=[T_ISO_VALUE]
)

# ==============================================================================
# 4. EMBEDDED BOUNDARY (EB) CYLINDER GEOMETRY
# ==============================================================================
dir_map = {0: (1.0, 0.0, 0.0), 1: (0.0, 1.0, 0.0), 2: (0.0, 0.0, 1.0)}
cylinder_dir = dir_map[EB_DIRECTION_AXIS]

burner_geometry = pv.Cylinder(
    center=EB_CENTER,
    direction=cylinder_dir,
    radius=EB_RADIUS,
    height=EB_HEIGHT,
    resolution=120,
)

# ==============================================================================
# 5. DESKTOP INTERACTIVE RENDERING
# ==============================================================================
print("[INFO] Launching PyVista Plotter...")
plotter = pv.Plotter()

# Render internal volume (T >= 800 K) with full colormap scalar range
plotter.add_mesh(
    clipped_volume,
    scalars="temp",
    cmap="inferno",
    clim=[T_AMBIENT, temp_max],  # Maps 300 K to min color, T_max to bright yellow
    smooth_shading=True,
    show_scalar_bar=True,
    scalar_bar_args={
        "title": "Temperature [K]",
        "vertical": True,
        "position_x": 0.85,
        "position_y": 0.1,
    },
)

# Render semi-transparent T = 800 K outer flame boundary shell
plotter.add_mesh(
    flame_shell,
    color="white",
    opacity=0.2,
    smooth_shading=True,
    label="Flame Boundary (T = 800 K)",
)

# Render Burner Nozzle Wall Mesh (EB Geometry)
plotter.add_mesh(
    burner_geometry,
    color="slategrey",
    metallic=0.7,
    roughness=0.3,
    show_edges=False,
    opacity=1.0,
    label="Burner Nozzle (EB)",
)

# Scene Configuration
plotter.set_background("#1e1e24")
plotter.add_axes(line_width=2)
plotter.add_bounding_box(color="white", line_width=0.8)
plotter.show_grid(color="gray", location="outer")

# View Setup
plotter.view_isometric()
plotter.reset_camera()

# Display Interactive Window
plotter.show(title="PeleLMeX LES - Internal Flame Volume (T >= 800 K)")