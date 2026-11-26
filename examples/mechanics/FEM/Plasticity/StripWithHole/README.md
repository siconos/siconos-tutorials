# Strip with a circular hole: Von Mises plasticity (plane stress)

Classical elastic–perfectly plastic test: a strip with a central circular hole pulled in
tension. By symmetry a quarter is modelled (`mesh_data/strip_with_circular_hole_72.msh`,
41 nodes, 60 T3 elements; `strip_with_circular_hole.geo` is the gmsh geometry):

- x ∈ [0, 0.5], y ∈ [0, 1], hole of radius 0.25 centred at the origin;
- bottom edge (`Dirichlet BC1`): v_y = 0; left edge (`Dirichlet BC2`): v_x = 0;
- top edge (`Prescribed displacement`): v_y = 100 (constant);
- E = 210e3, ν = 0.3, ρ = 7.8e-9, σ_y = 210, plane stress, lumped mass;
- h = 1e-6, θ = 1, 100 steps (about 0.4 s).

The reaction of the top edge grows linearly, then saturates (56.51 at the end, against the
net-section limit load σ_y · 0.25 = 52.5, overestimated by the coarse mesh); the plastic zone
starts at the root of the hole (step 4) and spreads across the ligament.

## Model

- `SolidLinearTIDS` (unknowns v and the element stresses, `MassType::lumped`,
  `Material(ρ, E, ν, AnalysisType2D::plane_stress)`);
- one `StressLinearTIR` per element, `C = elementStressSelection(e)` (its output is the
  physical stress of the element), with a `VonMisesPlasticityNSL(σ_y, plane_stress)`;
- `fem::integrators::MoreauJeanGOSI` and `fem::nonsmooth_formulations::Plasticity`: the global
  problem is condensed on the stresses and solved by the numerics `PlasticityProblem`
  (Von Mises, global semi-smooth Newton `PLASTICITY_2D_NSN`, tolerance 1e-10).

The scheme, the condensation and the validation are described in
`docs/devel/FEM_PLASTICITY_SCHEME.md` (siconos repository). The results agree with the
pytripopfem code (same mesh, native T3 backend, ipopt) to 4e-6 on the displacements and
4.5e-6 on the reaction.

The condensed matrix is ill-conditioned in this quasi-static regime: NSGS
(`PLASTICITY_2D_NSGS`) needs up to ~7e4 sweeps per step (about 2 min for the 100 steps),
the semi-smooth Newton 1 iteration per elastic step and 5–13 per plastic step.

## Outputs

- `strip_with_hole.dat`, one line per step: t, u_y of the top edge, reaction of the top
  edge, max σ_vm/σ_y, number of plastic elements, max norm of the plastic strain (compared
  with `strip_with_hole.ref`);
- `strip_with_hole_elements.dat`, final state per element: centroid (x, y), σxx, σyy, σxy,
  σ_vm, plastic strain (εxx, εyy, γxy);
- for the visualization: `strip_with_hole_nodes.dat` (x, y, dofs of u_x and u_y),
  `strip_with_hole_triangles.dat` (node indices) and the histories, one line per step from
  the initial state: `strip_with_hole_u.dat` (displacements, dof order),
  `strip_with_hole_sigma.dat`, `strip_with_hole_eps_p.dat` (physical stresses and plastic
  strains, 3 per element, Voigt order).

## Visualization (pyvista)

`strip_with_hole_viz.py` (adapted from `viz_pyvista.py` of pytripopfem) shows, on the
deformed mesh, the von Mises stress, the norm of the plastic strain and the displacement
magnitude, with colour scales fixed over the whole history. Run it in the directory of the
simulation outputs (or give it with `--dir`):

```bash
python3 strip_with_hole_viz.py                    # last step -> strip_with_hole_100.png
python3 strip_with_hole_viz.py --step 6 --warp 20 # first plastic steps, larger warp
python3 strip_with_hole_viz.py --gif              # all steps -> strip_with_hole.gif (imageio)
python3 strip_with_hole_viz.py --show             # interactive window
```

It needs `pyvista` (and `imageio` for `--gif`); it is not run by the tests.

Load–displacement curve:

```python
import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("strip_with_hole.dat")
plt.plot(d[:, 1], d[:, 2])
plt.axhline(210 * 0.25, ls="--", label="limit load")
plt.xlabel("u_top")
plt.ylabel("reaction")
plt.legend()
plt.show()
```
