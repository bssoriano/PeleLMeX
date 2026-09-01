#!/usr/bin/env python3
"""
3D flame shape, in the style of Fig. 5 of Mohit et al. (2026): the
INSTANTANEOUS T=800K flame boundary (an actual geometric iso-contour, via
marching cubes -- not a soft/blurred volume-rendered glow), with the
internal volume colored by temperature via a handful of further nested
iso-surfaces (800 -> ~2000 K) shown through the translucent outer shell.
Renders a single snapshot (like the paper's Fig. 5, which is explicitly
instantaneous, not time-averaged like flame_shape_history.py's Table
2/3-style metrics).

Uses PyVista/VTK, off-screen (no display/GPU needed -- verified working
headless on this machine via VTK's bundled software rasterizer).

Data pipeline (yt) is the same as flame_shape_history.py: full-domain
covering_grid from the plt file, cropped in-memory to a box around the
flame (found via the largest-connected-component of the T>=800K mask,
same rationale as that script -- this case's startup ignition pulse
otherwise leaves a disconnected hot blob in the field) -- the domain is
~63D long and the flame only occupies a few meters of it.

Axis convention for this case (see flame_shape_history.py for how this was
determined): y=downstream/crossflow, x=jet direction, z=spanwise.

Usage:
    python3 plot_flame_3d.py [plt_file]
"""
import sys
import numpy as np
import matplotlib.cm as cm
from scipy.ndimage import label
import pyvista as pv
import yt

yt.funcs.mylog.setLevel(50)
pv.OFF_SCREEN = True

PLT_FILE = "plt39127"      # instantaneous snapshot to render
T_FLAME = 800.0            # K, outer iso-contour (the flame boundary)
T_MAX = 2000.0             # K, colormap/innermost-layer upper bound
GRID_LEVEL = 2              # AMR level resampled onto (0=coarsest)
MARGIN = 0.15               # m, padding added around the flame's bounding box

# Nested iso-surfaces from the boundary inward, each more opaque than the
# last, so the translucent 800K shell reveals the hotter interior (found by
# trial and error -- see chat history -- that MANY thin layers with a
# smooth opacity ramp wash out to a grayish blend where the wrinkled
# surface causes lots of overlap; a handful of well-separated levels with
# a deliberately front-loaded opacity ramp reads much more like Fig. 5's
# solid-looking, clearly-shaded flame).
LAYER_LEVELS = [800, 1050, 1300, 1550, 1800, 1980]
LAYER_OPACITY = [0.45, 0.55, 0.65, 0.78, 0.9, 1.0]

OUT_PNG = "flame_3d_render.png"


def isolate_attached_flame(mask):
    """Keep only the largest connected T>=800K blob -- the real attached
    flame (tens of thousands of cells), as opposed to disconnected
    transient/noise blobs (far fewer cells). See flame_shape_history.py
    for the fuller rationale (this case's startup ignition pulse)."""
    labeled, n = label(mask)
    if n <= 1:
        return mask
    counts = np.bincount(labeled.ravel())
    counts[0] = 0
    return labeled == np.argmax(counts)


def load_cropped_temperature(plt_file):
    ds = yt.load(plt_file)
    level = min(GRID_LEVEL, ds.max_level)
    dims = (ds.domain_dimensions * ds.refine_by**level).astype(int)
    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)
    T = cg["temp"].v
    has_volfrac = ("boxlib", "volFrac") in ds.field_list
    vfrac = cg["volFrac"].v if has_volfrac else None

    le = ds.domain_left_edge.v
    re_ = ds.domain_right_edge.v
    dx = (re_ - le) / dims
    x = le[0] + (np.arange(dims[0]) + 0.5) * dx[0]
    y = le[1] + (np.arange(dims[1]) + 0.5) * dx[1]
    z = le[2] + (np.arange(dims[2]) + 0.5) * dx[2]

    mask = T >= T_FLAME
    if vfrac is not None:
        mask &= vfrac > 1e-3
    if not np.any(mask):
        raise RuntimeError(f"No cell in {plt_file} reaches T_FLAME={T_FLAME}K.")
    mask = isolate_attached_flame(mask)
    # zero out everything that isn't the attached flame so contour() at
    # T_FLAME can't pick up the disconnected transient/noise blobs
    T_clean = np.where(mask, T, 0.0)

    ix = np.where(np.any(mask, axis=(1, 2)))[0]
    iy = np.where(np.any(mask, axis=(0, 2)))[0]
    iz = np.where(np.any(mask, axis=(0, 1)))[0]
    x_lo, x_hi = x[ix[0]] - MARGIN, x[ix[-1]] + MARGIN
    y_lo, y_hi = y[iy[0]] - MARGIN, y[iy[-1]] + MARGIN
    z_lo, z_hi = z[iz[0]] - MARGIN, z[iz[-1]] + MARGIN
    print(f"Flame bounding box (+/-{MARGIN}m margin): "
          f"x=[{x_lo:.3f},{x_hi:.3f}] y=[{y_lo:.3f},{y_hi:.3f}] z=[{z_lo:.3f},{z_hi:.3f}]")

    icx = np.where((x >= x_lo) & (x <= x_hi))[0]
    icy = np.where((y >= y_lo) & (y <= y_hi))[0]
    icz = np.where((z >= z_lo) & (z <= z_hi))[0]
    T_crop = T_clean[np.ix_(icx, icy, icz)]
    print(f"Cropped field shape {T_crop.shape}, T range [{T_crop.min():.1f}, {T_crop.max():.1f}] K")

    spacing = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
    origin = (x[icx[0]], y[icy[0]], z[icz[0]])
    t_sim = ds.current_time.to_value("s")
    return T_crop, origin, spacing, t_sim


def main():
    plt_file = sys.argv[1] if len(sys.argv) > 1 else PLT_FILE
    T_crop, origin, spacing, t_sim = load_cropped_temperature(plt_file)

    grid = pv.ImageData(dimensions=T_crop.shape, spacing=spacing, origin=origin)
    grid.point_data["temperature"] = T_crop.flatten(order="F")

    pl = pv.Plotter(off_screen=True, window_size=[1400, 900])
    pl.set_background("white")

    cmap = cm.get_cmap("hot")
    for i, (lev, op) in enumerate(zip(LAYER_LEVELS, LAYER_OPACITY)):
        surf = grid.contour(isosurfaces=[lev], scalars="temperature")
        if surf.n_points == 0:
            continue
        rgba = cmap((lev - T_FLAME) / (T_MAX - T_FLAME))
        pl.add_mesh(
            surf, color=rgba[:3], opacity=op, smooth_shading=True,
            show_scalar_bar=False,
        )
    # a dummy zero-opacity actor just to get a temperature-labeled scalar
    # bar matching the layer colors, without letting it visually overlap
    # the (mesh-colored, not scalar-mapped) surfaces above
    outline = grid.outline()
    pl.add_mesh(
        outline, scalars=np.linspace(T_FLAME, T_MAX, outline.n_points),
        cmap="hot", opacity=0.0, show_scalar_bar=True,
        scalar_bar_args=dict(title="Temperature [K]", color="black"),
    )

    # NOTE: overlapping translucent surfaces need depth peeling for correct
    # blending -- without it, regions where several layers overlap (common
    # on a highly wrinkled turbulent surface) render as a muddy gray
    # instead of the correct blended color (see chat history).
    pl.enable_depth_peeling(number_of_peels=16, occlusion_ratio=0.0)

    pl.add_axes(color="black")
    # NOTE: `pl.camera_position = "iso"` (and any absolute camera.position/
    # focal_point assignment tried in its place) reproducibly makes meshes
    # render fully invisible here, for reasons not tracked down -- only
    # reset_camera() followed by RELATIVE azimuth/elevation adjustments was
    # reliable, so that's what's used below.
    pl.reset_camera()
    pl.camera.azimuth += 25
    pl.camera.elevation += 15

    pl.screenshot(OUT_PNG)
    print(f"Saved {OUT_PNG}  (t = {t_sim:.4f} s, from {plt_file})")


if __name__ == "__main__":
    main()
