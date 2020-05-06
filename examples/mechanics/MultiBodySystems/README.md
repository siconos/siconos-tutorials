Siconos mechanisms examples
---------------------------

[![A video of a watch escapment simulation](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/blob/master/examples/mechanics/MultiBodySystems/WatchEscapment/WatchEscapementMechanism.png)](https://www.youtube.com/watch?v=FhxgG4FdAgQ)


This directory contains examples of MultiBody Systems using OpenCascade for dealing with B-Rep (CAD step files). Contrary to the old interface in `siconos/mechanisms`, we use direclty hdf5 file to store the scene and the simulation results.

Model Definition and usage
--------------------------

### Watch escapement mechanism

The multibody system is defined in a python file, for instance `occ_watch_escapement.py`. To run the simulation, you have to

Compile the plugin that defines the spring forces:

	siconos --build-plugins

Run the python script:

	python occ_watch_escapement.py


The simulation results can be viewed from the hdf5 file :

    siconos_vview occ_watch_escapement.hdf5

