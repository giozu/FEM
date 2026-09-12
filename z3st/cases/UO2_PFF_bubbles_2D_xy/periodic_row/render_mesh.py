import sys

import pyvista as pv

vtk_path = sys.argv[1]
out_png = sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else ""

g = pv.read(vtk_path)

plotter = pv.Plotter(off_screen=True, window_size=[1100, 700])
plotter.add_mesh(g, color="lightsteelblue", show_edges=True, edge_color="gray", line_width=0.4)
plotter.view_xy()
plotter.camera.parallel_projection = True
bounds = g.bounds
cx = 0.5 * (bounds[0] + bounds[1])
cy = 0.5 * (bounds[2] + bounds[3])
half_h = 0.5 * (bounds[3] - bounds[2]) * 1.08
plotter.camera.SetFocalPoint(cx, cy, 0)
plotter.camera.SetPosition(cx, cy, half_h * 20)
plotter.camera.SetClippingRange(half_h * 0.1, half_h * 1000)
plotter.camera.SetParallelScale(half_h)
plotter.background_color = "white"
if title:
    plotter.add_text(title, position="upper_edge", font_size=12, color="black")
plotter.screenshot(out_png)
print(f"Wrote {out_png}  bounds={bounds}")
