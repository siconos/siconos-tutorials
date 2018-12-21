Siconos mechanisms examples
---------------------------

[![A video of a watch escapment simulation](./MultiBodySystems/WatchEscapment/WatchEscapementMechanism.png)](https://www.youtube.com/watch?v=FhxgG4FdAgQ)



This directory contains examples of MultiBody Systems using OpenCascade for dealing with B-Rep (CAD step files).

Model Definition and usage
--------------------------

The multibody system is defined in a python file, for instance occ_watch_escapement.py

    python occ_watch_escapement.py


The simulation results can also be viewed from the hdf5 file :

    siconos_vview occ_watch_escapement.hdf5
