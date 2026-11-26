# Strip footing on a Mohr–Coulomb soil (plane strain)

Bearing capacity of a strip footing on a cohesive-frictional, weightless soil. By symmetry
half of the problem is modelled (`mesh_data/half_footing_2x2_b02_size01.msh`, 356 nodes,
639 T3 elements, from pytripopfem):

- soil [0, 2] × [0, 2]; bottom edge (`Dirichlet BC1`): v_y = 0; lateral edges
  (`Dirichlet BC2`, `Dirichlet BC3`, the left one being the symmetry axis): v_x = 0;
- smooth rigid footing of half-width 0.2 on the top edge (`Forces`): v_y = −0.01 (constant,
  displacement control);
- E = 1e7 Pa, ν = 0.3, ρ = 1.5e3 kg/m³, plane strain, lumped mass;
- Mohr–Coulomb: cohesion c = 2 kPa, friction angle φ = 25°, dilatancy angle ψ = φ
  (associated flow) by default;
- h = 0.05, θ = 1, 40 steps (about 13 s).

The footing reaction grows and saturates near −11 900 (half footing). Prandtl's bearing
capacity of a weightless soil, q_u = c N_c (N_c = 20.72 for φ = 25°), gives −8 288 for the
half footing: the constant-strain T3 elements (stiff in plane strain plasticity) and the small
domain overestimate it. The plastic zone starts under the footing (step 1) and spreads to the
whole domain, with a heave next to the footing.

## Model

- `SolidLinearTIDS` (unknowns v and the element stresses, `MassType::lumped`,
  `Material(ρ, E, ν, AnalysisType2D::plane_strain)`);
- one `StressLinearTIR` per element, `C = elementStressSelection(e)` (its output is the
  physical stress of the element), with a
  `DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain(c, φ, ψ)`: R ≤ −p sin φ + c cos φ (R the
  radius of the Mohr circle, p the in-plane mean stress), i.e. the cone ‖σ_T‖ ≤ η σ_N + c' of
  the numerics `PlasticityProblem` with η = sin φ, c' = √2 c cos φ, and the dilatancy parameter
  θ = sin φ − sin ψ;
- `fem::integrators::MoreauJeanGOSI` and `fem::nonsmooth_formulations::Plasticity`: the global
  problem is condensed on the stresses and solved by the global semi-smooth Newton
  `PLASTICITY_2D_NSN` (tolerance 1e-10, 7–13 iterations per step; NSGS reaches an error of
  3e-4 after 2e4 sweeps on these problems).

The results agree with the pytripopfem code (same mesh, native T3 backend, its semi-smooth
Newton solver) to round-off (relative differences ≤ 7e-12 on the displacements, 3e-11 on
the stresses, 4e-13 on the footing reaction), for ψ = φ and for ψ = 10° (non-associated
flow), see `docs/devel/FEM_PLASTICITY_SCHEME.md`
(siconos repository).

Usage: `./footing [mesh [ψ in degrees]]`, e.g. `./footing
mesh_data/half_footing_2x2_b02_size01.msh 10` (the reference comparison is done for the
default case only).

## Outputs

- `footing.dat`, one line per step: t, u_y of the footing, footing reaction, max of the
  yield function f / (c cos φ) (f = R + p sin φ − c cos φ ≤ 0), number of plastic elements,
  max norm of the plastic strain (compared with `footing.ref`);
- for the visualization: `footing_nodes.dat` (x, y, dofs of u_x and u_y),
  `footing_triangles.dat` (node indices) and the histories, one line per step from the
  initial state: `footing_u.dat` (displacements, dof order), `footing_sigma.dat`,
  `footing_eps_p.dat` (physical stresses and plastic strains, 3 per element, Voigt order).

## Visualization (pyvista)

`footing_viz.py` shows, on the deformed mesh, the yield function f / (c cos φ) (−1 for a
stress-free element, 0 on the yield surface), the norm of the plastic strain and the
displacement magnitude, with colour scales fixed over the whole history. Run it in the
directory of the simulation outputs (or give it with `--dir`):

```bash
python3 footing_viz.py                    # last step -> footing_040.png
python3 footing_viz.py --step 2           # onset of plasticity under the footing
python3 footing_viz.py --gif              # all steps -> footing.gif (imageio)
python3 footing_viz.py --show             # interactive window
```

`--cohesion` and `--phi` must match the simulation. It needs `pyvista` (and `imageio` for
`--gif`); it is not run by the tests.

Load–displacement curve:

```python
import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("footing.dat")
plt.plot(-d[:, 1], -d[:, 2])
plt.axhline(8288, ls="--", label="Prandtl (half footing)")
plt.xlabel("settlement")
plt.ylabel("footing reaction")
plt.legend()
plt.show()
```
