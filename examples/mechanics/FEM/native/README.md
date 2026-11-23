# Tests of finit element model implemented in Siconos

Available examples :

- T3...*.cpp : 2D mesh (square domain) with T3 elements.
- TH4.cpp : 3D (box) mesh with TH4 elements
- Hertz_... : 2D disk (T3 elems), placed on a plane and then compressed in one direction.

Some utilities functions are available in src/native_fem_utils.*


All examples are built using gmsh input files.

Results and post-processing :

```
siconos XX.cpp
````

results in:
- XX.dat to be compared with some reference file (check in cpp file for the name of the ref file)
- outputs/YY.py (mesh file) if YY.msh has been used
- outputs/XX_displacement.py to be used as input for meshio processing:

```
python meshio_prepost.py --mesh_file="outputs/mesh_data/YY.msh" --simulation_name="XX"
```
will generate vtk files in ./vtk/XX... than can be viewed with Paraview.
