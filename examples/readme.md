Siconos examples
================

This directory contains Siconos examples programs, sorted by application area.

All examples can be executed independently.

They require a proper install of Siconos (see [Download and install Siconos](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/install_guide/index.html).


Run a simulation
----------------

Try for example :

```shell
cd mechanics/BouncingBall
siconos BouncingBallTS.cpp

siconos BouncingBallTS.py
```


Build and run all examples
--------------------------


```shell
mkdir build
cd build
cmake <PATH_TO_SICONOS_TUTORIAL>/examples -Dsiconos_DIR=<PATH_TO_SICONOS_INSTALL>/share/siconos/cmake/
make -j <N>
make test -j <N>
```

with :

* PATH_TO_SICONOS_TUTORIAL : the path to your siconos-tutorial git repository
* PATH_TO_SICONOS_INSTALL : where you have installed Siconos software
* N : number of processes to be used for compilation.

You will get a report with the list of working/non working examples.

For an uptodate dashboard of continious integration process of Siconos examples, please check [Siconos cdash](https://cdash-tripop.inrialpes.fr/index.php?project=Siconos), build name "Siconos examples".
