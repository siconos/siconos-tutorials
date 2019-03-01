# Title

  This example uses old inpout files of LMGC77 created by V. Acary in his PhD Thesis (2001)

# How to

-  To parse the LMGC77 input file, run

       python parser_lmgc77_2D.py

-  To build 3D brick form 2D Finite element mesh run

     python make_3d_brick.py

- to run the simulation

      python arches.py

- to postprocess

     siconos_vexport --global-filter arches.hdf5

     siconos_vview --normalcone-ratio=0.01 --cf-scale=0.2 arches.hdf5 