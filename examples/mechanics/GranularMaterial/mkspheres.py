#!/usr/bin/env python

# generation of some spheres in a box
# with radii between radius_min and radius_max
# generation in a 20cm x 20cm 20cm box
# ./mkspheres.py <nb spheres>
# output :
#  - coordinates: coors-<nb generated spheres>.txt
#  - radii      : radii-<nb generated spheres>.txt


import sys

import numpy

if __name__ == "__main__":
    try:
        from pylmgc90 import pre

    except ImportError:
        print("pylmgc90 is not found !!!, exit")
        sys.exit(0)

radius_min = 0.001
radius_max = 0.0025
lx = 0.02
ly = 0.02
lz = 0.02

print(sys.argv)


if __name__ == "__main__":
    try:
        from pylmgc90 import pre
    except ImportError:
        print("pylmgc90 is not found !!!, exit")
        sys.exit(0)
    if len(sys.argv) > 1:
        nbp = int(sys.argv[1])
    else:
        nbp = 20
    radii = pre.granulo_Random(nbp, radius_min, radius_max)
    print(radii)
    print(pre.depositInBox3D(radii, lx, ly, lz))
    nbpl, coors, radii = pre.depositInBox3D(radii, lx, ly, lz)

    numpy.savetxt(f"coors-{nbpl}.txt", coors[0 : nbpl * 3])
    numpy.savetxt(f"radii-{nbpl}.txt", radii[0:nbpl])
