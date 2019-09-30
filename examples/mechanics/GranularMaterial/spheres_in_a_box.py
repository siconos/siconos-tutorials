#!/usr/bin/env python

# spheres in a box of size:
# mkspheres.lx*mkspheres.ly*mspheres.lz
# for n spheres:
# ./mkspheres.py <n>
# ./spheres_in_a_box.py coors-<n>.txt radii-<n>.txt

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.io.FrictionContactTrace import FrictionContactTraceParams
import siconos.numerics as Numerics
from siconos.mechanics.collision.bullet import SiconosBulletOptions
from math import pi
import numpy
import sys
import mkspheres
options = SiconosBulletOptions()
options.worldScale = 1000
options.contactBreakingThreshold = 0.0002
options.perturbationIterations = 0
options.minimumPointsPerturbationThreshold = 0
hstep = 0.001
theta = 0.50001
itermax = 1000
tolerance = 1e-7
lx = mkspheres.lx
ly = mkspheres.ly
lz = mkspheres.lz
margin_ratio=1.e-5
wthick = mkspheres.lz/10
zoffset = -wthick/2

margin_max = margin_ratio*mkspheres.radius_max

coors_filename=sys.argv[1]
radii_filename=sys.argv[2]

coors = numpy.loadtxt(coors_filename)
radii = numpy.loadtxt(radii_filename)

if len(radii.shape) == 0:
    radii=[radii]

nb_laid_particles=len(radii)

print(nb_laid_particles)

solver=Numerics.SICONOS_FRICTION_3D_NSGS
fileName = "spheres-in-a-box-{0}".format(nb_laid_particles)
title = "SpheresBox"
description = """
Spheres in a box, generation with lmgc90 granulo_Random
number of spheres: {0}
radius min       : {1}
radius max       : {2}
box size x       : {3}
box size y       : {4}
box size z       : {5}
Moreau TimeStepping: h={6}, theta={7}
One Step non smooth problem: {8}, maxiter={9}, tol={10}
""".format(nb_laid_particles,
           mkspheres.radius_min,
           mkspheres.radius_max,
           lx, ly, lz,
           hstep, theta, Numerics.solver_options_id_to_name(solver),
           itermax,
           tolerance)

mathInfo = ""

friction_contact_trace_params = FrictionContactTraceParams(
    dump_itermax=10000, dump_probability=None,
    fileName=fileName, title=title,
    description=description, mathInfo=mathInfo)

with MechanicsHdf5Runner(io_filename='siab-{0}.hdf5'.format(nb_laid_particles)) as io:
    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (lx, ly, wthick),
                           insideMargin=margin_max, outsideMargin=margin_max)
    io.add_primitive_shape('Wall', 'Box', (lx, ly, wthick),
                           insideMargin=margin_max, outsideMargin=margin_max)

    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, zoffset])
    io.add_object('wall1', [Contactor('Wall')],
                  translation=[-lx/2-wthick/2, 0, lz/2+zoffset],
                  orientation=([0, 1, 0], pi/2.))
    io.add_object('wall2', [Contactor('Wall')],
                  translation=[lx/2+wthick/2, 0, lz/2+zoffset],
                  orientation=([0, 1, 0], pi/2.))
    io.add_object('wall3', [Contactor('Wall')],
                  translation=[0, -ly/2-wthick/2, lz/2+zoffset],
                  orientation=([1, 0, 0], pi/2.))
    io.add_object('wall4', [Contactor('Wall')],
                  translation=[0, ly/2+wthick/2, lz/2+zoffset],
                  orientation=([1, 0, 0], pi/2.))

    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.)
    for i in range(nb_laid_particles):
        rad=radii[i]
        margin=rad*margin_ratio
        io.add_primitive_shape('Sphere-{0}'.format(i), 'Sphere', (rad,),
                               insideMargin=margin, outsideMargin=margin)
        mass = (4./3)*pi*(radii[i]**3)*2320
        I_ = (2./5)*mass*(radii[i]*radii[i])
        io.add_object('sphere-{0}'.format(i), [Contactor('Sphere-{0}'.format(i))],
                      translation=[coors[3*i], coors[3*i+1], coors[3*i+2]],
                      velocity=[0, 0, 0, 0, 0, 0],
                      inertia=[I_, I_, I_],
                      mass=mass)

with MechanicsHdf5Runner(io_filename='siab-{0}.hdf5'.format(nb_laid_particles), mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=True,
           options=options,
           face_class=None,
           edge_class=None,
           t0=0,
           T=8.00,
           h=hstep,
           multipoints_iterations=False,
           theta=theta,
           Newton_max_iter=1,
           set_external_forces=None,
           friction_contact_trace=True,
           friction_contact_trace_params=friction_contact_trace_params,
           constraint_activation_threshold=0.00001,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=itermax,
           tolerance=tolerance,
           numerics_verbose=False,
           explode_Newton_solve=True,
           output_frequency=None)
