# siconos-tutorials 2026.11 Release Notes

**Release Date:** November 2026
**Siconos version:** 2026.11 (use the same version of siconos and siconos-tutorials)
**Previous Release:** 4.4.0 (September 2021)
**Commits:** about 270 since 4.4.0

---

## Overview

This release follows the modernization of Siconos between 4.4.0 and 2026.11: C++20 and the
`siconos::` namespaces, Eigen-based algebra, pybind11 Python bindings (SWIG removed), user
functions given as C++ lambdas or Python callables instead of plugins loaded by name, and the
renaming of the Lagrangian and Newton-Euler relations. All the C++ and Python examples run by
the tests were ported to this API, their references checked, and the tests made independent of
the platform (macOS arm64, Linux x86_64).

The notes of Siconos itself are in `docs/devel/RELEASE_2026.11.md` of the siconos repository.

---

## Porting to the Siconos 2026.11 API

- **C++**: namespaces (`siconos::modeling`, `siconos::simulation`, ...), `std::shared_ptr`
  instead of `SP::`, `SiconosMatrix`/`SiconosVector` (Eigen) instead of `SimpleMatrix`,
  matrices initialized explicitly (`setZero()`, Eigen does not initialize).
- **Python**: pybind11 modules (`siconos.modeling`, `siconos.simulation`,
  `siconos.integrators`, `siconos.nonsmooth_formulations`, `siconos.numerics`, ...),
  `copy_t`/`alias_t` arguments of the constructors and setters, numpy arrays in column-major
  order for matrices.
- **User functions**: lambdas (`setCompute...Function`) instead of plugins loaded by name;
  `cppimport` removed (addons built by the `siconos` script).
- **Relations**: descriptive names (`contactMatrix`, `constraintJacobian`,
  `Lagrangian2dVelocityR`, `NewtonEulerVelocityR`, ...), new `computeh` API of `NewtonEulerR`.
- **Numerics**: new solver names and options header (`FrictionContact_options.h`).
- **mechanics_run**: options given by a `run_options` dictionary (former keyword arguments
  and backend `*_class` options removed); `explode_computeOneStep_in_python` and
  `explode_computeOneStepNSProblem_in_python` options; energy and work outputs; restart.

---

## New Examples

- **Rock fall** (`mechanics/RockFall`): rock fall from a digital elevation model (raster file,
  with the Frémond law), on a height map, on a forest, rock fall propagation with restart,
  rolling example, chute with a choice of rock shape.
- **FEM** (`mechanics/FEM`): native T3 and TH4 elements, T3 with collision detection (with and
  without friction), Hertz contact, T3 with `MoreauJeanGOSI`, Euler-Bernoulli beam and ball
  against a beam, impacting bar with sparse matrices (C++ and Python), plastic footing
  (`FEM/Plasticity/T3-Footing`).
- **Frémond impact law**: `BouncingBallTS-FremondNSL`, `box_stacks_FremondNSL.py`,
  `rockfall_from_dem_FremondNSL.py`, rock protection wall with energy outputs.
- **Boundary conditions**: `mechanics/BallOnMovingPlane` (C++ and Python), reactions of the
  boundary conditions in the outputs.
- **Schemes**: `RockingBlock` with both projection updates (incremental and `from_origin`),
  `SliderCrankD1MinusLinear` with both impact rules, stiff linear system for Lsodar,
  `ColumnOfBeadsTS-MoreauJeanGOSI`, pendulum with Lsodar and the alpha scheme.
- **vKernel** (built when siconos is installed with `WITH_vkernel=ON`):
  `BouncingBallTS-VkernelMoreauJean`, `BouncingBallTS-VkernelFremondNSL`.
- Others: transport cable (`mechanics/TransportCable`), friction oscillator in Python,
  multicylinder and relative orientation tests, primitive soup in multiple runs, removal of
  static contactors, granular chute with vibrators.

---

## Tests

- One ctest test per example: `build_<dir>_<example>` and `run_<dir>_<example>` for C++,
  `py_<dir>_<script>` for Python; the tests use the components found in the siconos
  installation.
- Short runs: the tests set `SICONOS_EXAMPLES_SHORT_RUN`, the longest examples read it to
  simulate a shorter time; the full simulation is run otherwise.
- References: comparisons relative to the size of the data, groups of columns with their own
  tolerance when some quantities are amplified (accelerations and forces of the
  generalized-alpha scheme, event-driven examples), checks of properties instead of
  references for the chaotic cases (`BulletBouncingBox`).
- Ubuntu 26.04 (GCC 15, Python 3.14, numpy 2.5) fixes; formatting and lint with pre-commit
  hooks (clang-format, ruff, shellcheck); CI on the `main` branch.
- Status with Siconos 2026.11 on macOS arm64: the 332 tests pass. On Ubuntu 26.04 x86_64 they
  passed before the last fixes; under x86_64 emulation, `mechanics/FEM/native/T3-GOSI` may
  time out.

---

## Removed

- `cppimport` (replaced by the addons of the `siconos` script).
- `mechanics/NewtonEuler/BouncingBallNETS.py`.
- FEM plasticity examples `plasticBernoulliBeam` and `ballAgainstPlasticBernoulliBeam`
  (`mechanics/FEM/native`), removed for this release while the FEM plasticity formulation of
  Siconos is reviewed.

---

## Known Issues

- **Python control examples are not updated** (`examples/control/*/*.py`): they use the
  former `siconos.kernel`/`siconos.control` modules, former constructors and numerics
  bindings that do not exist (`VI`, `MixedComplementarityProblem2`). They are not run by the
  tests: `examples/CMakeLists.txt` scans only the top of `examples/control`, so none of the
  control examples is a test. The C++ control examples build and run, but some references
  are old: small round-off differences for the observers and `SMC*`, and
  `TwoLinkManipulator_{Flexible, MultiConstraints, MultiConstraintsFlexible}` differ from
  their 2018 references (10⁻⁷ at the first steps, O(1) after the impacts).
- **Jupyter notebooks are not updated**:
  - `siconos-notebooks/st01_diode_bridge`, `st02_bouncing_ball`, `st03_bouncing_ball_ed`
    call constructors without the `copy_t`/`alias_t` argument;
  - `siconos-notebooks/advanced_01_impacting_bar` and
    `examples/electronics/DiodeBridge/DiodeBridge.ipynb` import `siconos.kernel`;
  - `siconos-notebooks/advanced_02_guitar_string` needs `numpywrappers` from
    `examples/mechanics/Music`;
  - `ControlExamples`, `ElecExamples` and `MechanicsExamples` plot the results of the
    examples, which must be run first.
  The Python scripts of the same examples are up to date.
- `Mechanisms/*` examples cannot run: the Python bindings of the siconos `mechanisms`
  component are not implemented.
- OCC examples `occ_watch_escapement` and `occ_slider_crank` use C plugins loaded by name,
  removed from Siconos.
- RockFall scripts need `trimesh`; `simu_rockfall_biplan` ignores the short-run mode.

---

## Contributors

Vincent Acary, Franck Pérignon, Olivier Goury, Samuel Heidmann, Maurice Brémond, F. Bourrier
(rock fall rolling example). See [AUTHORS](AUTHORS) and [credits.md](credits.md).

---

## Download

```bash
git clone https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials.git
cd siconos-tutorials
git checkout 2026.11

# with siconos 2026.11 installed
cmake -S examples -B build
ctest --test-dir build
```

**License**: Apache License 2.0
**Copyright**: INRIA
