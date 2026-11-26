# Tests of finite element models implemented in Siconos

Examples are grouped by model, one directory each:

- `T3_square/` : 2D mesh (square domain) with T3 elements (`T3*.cpp`).
  Utility functions are in `T3_square/src/native_fem_utils.*`
  (built automatically by the `siconos` script).
- `Hertz/` : 2D disk (T3 elements), placed on a plane and then compressed in one direction.
- `TH4/` : 3D (box) mesh with TH4 elements.
- `BernoulliBeam/` : elastic Euler-Bernoulli beam, alone and impacted by a ball.
  `BeamFEMForceField.scn.xml` is the equivalent SOFA scene, for comparison.
- `PlasticBernoulliBeam/` : elasto-plastic Euler-Bernoulli beam, alone and impacted by a ball.

The 2D/3D examples are built from gmsh input files, stored in the `mesh_data/`
directory of each example.

Post-processing scripts (`meshio_prepost.py`, `meshio_prepost_beam.py`,
`matplotlib_post_2d.py`) are shared and stay in this directory.

## Results and post-processing

From an example directory:

```
siconos XX.cpp
```

results in:
- XX.dat to be compared with some reference file (check in cpp file for the name of the ref file)
- outputs/YY.py (mesh file) if YY.msh has been used
- outputs/XX_displacement.py to be used as input for meshio processing:

```
python ../meshio_prepost.py --mesh_file="outputs/mesh_data/YY.msh" --simulation_name="XX"
```
will generate vtk files in ./vtk/XX... that can be viewed with Paraview.
