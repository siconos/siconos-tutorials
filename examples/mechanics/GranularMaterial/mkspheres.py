#!/usr/bin/env python

# generation of some spheres in a box
# with radii between radius_min and radius_max
# generation in a 20cm x 20cm 20cm box
# ./mkspheres.py <nb spheres>
# output :
#  - coordinates: coors-<nb generated spheres>.txt
#  - radii      : radii-<nb generated spheres>.txt




import numpy
from pylmgc90 import pre
import sys

radius_min = 0.001
radius_max = 0.0025
lx = .2
ly = .2
lz = .2

if __name__ == "__main__":
    nbp = int(sys.argv[1])
    
    radii = pre.granulo_Random(nbp, radius_min, radius_max)
    [nbpl, coors] = pre.depositInBox3D(radii, lx, ly, lz)

    numpy.savetxt('coors-{0}.txt'.format(nbpl), coors[0:nbpl*3])
    numpy.savetxt('radii-{0}.txt'.format(nbpl), radii[0:nbpl])
