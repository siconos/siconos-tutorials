"""Visualise the results of footing.cpp with pyvista.

Reads the outputs of the simulation (in --dir, the directory where footing was run): the
mesh (footing_nodes.dat, footing_triangles.dat) and the histories (footing_u.dat,
footing_sigma.dat, footing_eps_p.dat, footing.dat for the time), and shows, on the deformed
mesh, the Mohr-Coulomb yield function f / (c cos phi) (f = R + p sin(phi) - c cos(phi):
-1 for a stress-free element, 0 on the yield surface), the tensor norm of the plastic strain and
the displacement magnitude.

Adapted from viz_pyvista.py of pytripopfem (examples/Strip_with_hole).

Examples (without a display, images are saved):
    python3 footing_viz.py                    # last step -> png
    python3 footing_viz.py --step 10 --warp 20
    python3 footing_viz.py --gif              # all steps -> gif (needs imageio)
    python3 footing_viz.py --show             # interactive window (needs a display)
"""

import argparse
import os

import numpy as np
import pyvista as pv

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("--dir", default=".", help="directory of the simulation outputs")
parser.add_argument("--step", type=int, default=-1, help="step to plot (default: last)")
parser.add_argument(
    "--warp", type=float, default=10.0, help="displacement magnification"
)
parser.add_argument("--cohesion", type=float, default=2e3, help="cohesion c")
parser.add_argument("--phi", type=float, default=25.0, help="friction angle (degrees)")
parser.add_argument(
    "--eps-max",
    type=float,
    default=None,
    help="upper bound of the plastic strain colour scale (default: 98th percentile of the "
    "last step; the singular element at the edge of the footing would hide the rest)",
)
parser.add_argument(
    "--eps-decades",
    type=float,
    default=3.0,
    help="the plastic strain is shown with a log scale over this number of decades",
)
parser.add_argument(
    "--gif", action="store_true", help="write an animation of all steps"
)
parser.add_argument("--show", action="store_true", help="open an interactive window")
args = parser.parse_args()

prefix = os.path.join(args.dir, "footing")
if not os.path.exists(prefix + "_u.dat"):
    raise SystemExit(f"no outputs in {args.dir}: run footing first")
phi = np.radians(args.phi)
c_cos_phi = args.cohesion * np.cos(phi)

nodes = np.loadtxt(prefix + "_nodes.dat", ndmin=2)  # x, y, dof x, dof y
triangles = np.loadtxt(prefix + "_triangles.dat", dtype=int, ndmin=2)
u_hist = np.loadtxt(prefix + "_u.dat", ndmin=2)  # steps x ndof
sigma_hist = np.loadtxt(prefix + "_sigma.dat", ndmin=2)  # steps x 3 nel, Mandel
eps_hist = np.loadtxt(prefix + "_eps_p.dat", ndmin=2)  # steps x 3 nel, Mandel
t = (
    np.loadtxt(prefix + ".dat", ndmin=2)[:, 0]
    if os.path.exists(prefix + ".dat")
    else None
)
n_steps = u_hist.shape[0]
dof = nodes[:, 2:4].astype(int)

points = np.column_stack([nodes[:, 0], nodes[:, 1], np.zeros(len(nodes))])
cells = np.column_stack([np.full(len(triangles), 3), triangles]).ravel()
grid = pv.UnstructuredGrid(cells, np.full(len(triangles), pv.CellType.TRIANGLE), points)


def tensor_norm(eps):
    """norm of the strain tensors given in Mandel notation (exx, eyy, sqrt2 exy), one per
    element: the Euclidean norm, sqrt(exx^2 + eyy^2 + 2 exy^2)"""
    return np.linalg.norm(eps.reshape(-1, 3), axis=1)


def load(k):
    """mesh of step k with its fields, deformed by warp * u"""
    mesh = grid.copy()
    u = u_hist[k][dof]  # (n_nodes, 2)
    mesh.point_data["u"] = np.column_stack([u, np.zeros(len(u))])
    mesh.point_data["u_norm"] = np.linalg.norm(u, axis=1)
    s = sigma_hist[k].reshape(-1, 3)
    p = 0.5 * (s[:, 0] + s[:, 1])
    radius = np.hypot(0.5 * (s[:, 0] - s[:, 1]), s[:, 2] / np.sqrt(2.0))  # Mandel
    mesh.cell_data["yield_function"] = (
        radius + p * np.sin(phi) - c_cos_phi
    ) / c_cos_phi
    mesh.cell_data["eps_p_norm"] = tensor_norm(eps_hist[k])
    # log scale: values below the scale (elastic elements) are shown at its lower bound
    mesh.cell_data["eps_p_log"] = np.maximum(mesh.cell_data["eps_p_norm"], eps_floor)
    return mesh.warp_by_vector("u", factor=args.warp)


# colour ranges fixed over the whole history so that an animation is readable
eps_last = tensor_norm(eps_hist[-1])
eps_top = (
    args.eps_max if args.eps_max else max(float(np.percentile(eps_last, 98)), 1e-12)
)
eps_floor = eps_top * 10.0 ** (-args.eps_decades)
last = load(n_steps - 1)
panels = [
    ("yield_function", "Mohr-Coulomb f / (c cos phi)", "viridis", [-1.0, 0.0], "%.1f"),
    (
        "eps_p_log",
        "plastic strain tensor norm (log scale)",
        "magma",
        [eps_floor, eps_top],
        "%.0e",
    ),
    (
        "u_norm",
        "displacement norm",
        "plasma",
        [0.0, max(float(last["u_norm"].max()), 1e-12)],
        "%.1e",
    ),
]

plotter = pv.Plotter(
    shape=(1, len(panels)), off_screen=not args.show, window_size=(1500, 650)
)
for j, (_, title, *_) in enumerate(panels):
    plotter.subplot(0, j)
    plotter.add_text(title, font_size=10)
plotter.link_views()


def update(k):
    """draw step k in every panel (the meshes are replaced, camera kept)"""
    mesh = load(k)
    for j, (field, _, cmap, clim, fmt) in enumerate(panels):
        plotter.subplot(0, j)
        plotter.add_mesh(
            mesh,
            scalars=field,
            cmap=cmap,
            clim=clim,
            log_scale=(field == "eps_p_log"),
            show_edges=True,
            edge_color="grey",
            name="field",
            copy_mesh=True,
            scalar_bar_args={
                "title": field,
                "vertical": False,
                "fmt": fmt,
                "n_labels": 3,
            },
        )
    plotter.subplot(0, 0)
    time = f"  t = {t[k]:.3e}" if t is not None and k < len(t) else ""
    plotter.add_text(
        f"step {k}/{n_steps - 1}{time}   warp x{args.warp:g}",
        position="lower_left",
        font_size=9,
        name="label",
    )


update(n_steps - 1)  # first draw with the largest deformation, then fix the camera
for j in range(len(panels)):
    plotter.subplot(0, j)
    plotter.view_xy()

if args.gif:
    try:
        import imageio  # noqa: F401  (needed by pyvista for gif)

        out = "footing.gif"
        plotter.open_gif(out, fps=10)
        for k in range(n_steps):
            update(k)
            plotter.write_frame()
        print("written", out)
    except ModuleNotFoundError:
        out = "footing_frames"
        os.makedirs(out, exist_ok=True)
        for k in range(n_steps):
            update(k)
            plotter.screenshot(os.path.join(out, f"{k:03d}.png"))
        print(
            f"imageio not installed (pip install imageio): PNG frames written in {out}/"
        )
    plotter.close()
else:
    k = args.step % n_steps
    update(k)
    if args.show:
        plotter.show()
    else:
        out = f"footing_{k:03d}.png"
        plotter.screenshot(out)
        print("written", out)
