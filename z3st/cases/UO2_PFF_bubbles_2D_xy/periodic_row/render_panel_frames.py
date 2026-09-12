import os

import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6

# window fixed in absolute (x,y) for every panel
X0, X1 = -5e-6, 50e-6
Y0, Y1 = -2e-6, 32e-6
CX, CY = 0.5 * (X0 + X1), 0.5 * (Y0 + Y1)
HALF_H = 0.5 * (Y1 - Y0)
ASPECT = (X1 - X0) / (Y1 - Y0)

nucleation_steps = {200: 154, 100: 140, 60: 121, 45: 106, 36: 91, 33: 84, 30: 75}
d09_steps = {200: 156, 100: 142, 60: 123, 45: 108, 36: 93, 33: 86, 30: 77}

W = 900
H = int(W / ASPECT)

for Lx_um in nucleation_steps:
    for tag, step in [("nucleation", nucleation_steps[Lx_um]), ("d09", d09_steps[Lx_um]), ("final", 200)]:
        vtu = os.path.join(CASE, f"Lx_{Lx_um}", "output", f"simulation_{step:04d}.vtu")
        g = pv.read(vtu)
        plotter = pv.Plotter(off_screen=True, window_size=[W, H])
        plotter.add_mesh(g, scalars="Damage", cmap="turbo", clim=[0, 1],
                          show_scalar_bar=False, show_edges=False)
        plotter.view_xy()
        plotter.camera.parallel_projection = True
        plotter.camera.SetFocalPoint(CX, CY, 0)
        plotter.camera.SetPosition(CX, CY, HALF_H * 20)
        plotter.camera.SetClippingRange(HALF_H * 0.1, HALF_H * 1000)
        plotter.camera.SetParallelScale(HALF_H)
        plotter.background_color = "white"
        out = os.path.join(CASE, f"_panel_Lx{Lx_um}_{tag}.png")
        plotter.screenshot(out)
        plotter.close()
        print(f"wrote {out}")
