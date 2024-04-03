# A rock protection masonry wall againt rock fall hazards
   file:

	rock_protection_wall.py

![Protection masonry wall against rock fall](rock_protection_wall.jpg)

## various configurations inside the scripts
  - pyramid wall
  - wide wall
  - wide wall with buttresses
  - tall wall
  - tall wall with buttresses

## post-processing:

	siconos_vview rock_protection_wall.hdf5
or for paraview

	siconos_vexport --global-filter rock_protection_wall.hdf5

or for paraview in parrallel

	siconos_vexport --global-filter --gen-para-script=6  rock_protection_wall.hdf5 > par.sh; sh par.sh
