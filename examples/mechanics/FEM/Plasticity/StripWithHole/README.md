# Strip with a circular hole: Von Mises plasticity (plane stress)

Classical elastic–perfectly plastic test: a strip with a central circular hole pulled in
tension. By symmetry a quarter is modelled (`mesh_data/strip_with_circular_hole_72.msh`,
41 nodes, 60 T3 elements; `strip_with_circular_hole.geo` is the gmsh geometry):

- x ∈ [0, 0.5], y ∈ [0, 1], hole of radius 0.25 centred at the origin;
- bottom edge (`Dirichlet BC1`): v_y = 0; left edge (`Dirichlet BC2`): v_x = 0;
- top edge (`Prescribed displacement`): v_y = 100 (constant);
- E = 210e3, ν = 0.3, ρ = 7.8e-9, σ_y = 210, plane stress, lumped mass;
- h = 1e-6, θ = 1, 100 steps (20 in the tests: `SICONOS_EXAMPLES_SHORT_RUN`).

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
  (Von Mises, NSGS, tolerance 1e-8).

The scheme, the condensation and the validation are described in
`docs/devel/FEM_PLASTICITY_SCHEME.md` (siconos repository). The results agree with the
pytripopfem code (same mesh, native T3 backend, ipopt) to 3e-5 on the displacements and
7e-6 on the reaction.

The condensed matrix is ill-conditioned in this quasi-static regime: NSGS needs up to
~7e4 sweeps per step (about 25 s for the 20 steps of the test, 2 min for the full run).

## Outputs

- `strip_with_hole.dat`, one line per step: t, u_y of the top edge, reaction of the top
  edge, max σ_vm/σ_y, number of plastic elements, max norm of the plastic strain;
- `strip_with_hole_elements.dat`, final state per element: centroid (x, y), σxx, σyy, σxy,
  σ_vm, plastic strain (εxx, εyy, γxy).

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

Plastic zone:

```python
e = np.loadtxt("strip_with_hole_elements.dat")
plt.tricontourf(e[:, 0], e[:, 1], e[:, 5] / 210)
plt.colorbar(label="sigma_vm / sigma_y")
plt.axis("equal")
plt.show()
```
