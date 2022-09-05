#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk

import numpy as np
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cylinder
    R = 0.1
    L = 2.0
    io.add_primitive_shape('Cyl', 'Cylinder', (R, L))
    io.add_primitive_shape('Stick', 'Box', (R, L, R))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (20, 10, 1.0))
    
    # Definition of the ground shape
    io.add_primitive_shape('SmallBox', 'Box', (.1, .1, .1))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.0)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    mass_test = 1.0
    inertia_test = np.eye(3)

    inertia_test[0, 0] = 0.25*mass_test*R*R + 1/3.0*mass_test*L*L
    inertia_test[1, 1] = 0.5*mass_test*R*R
    inertia_test[2, 2] = 0.25*mass_test*R*R + 1/3.0*mass_test*L*L
    print(inertia_test)
    orientation_test = [0.14, 0.7, 0.7, 0]
    # io.add_object('cyl1', [Contactor('Cyl')],
    #               translation=[-2, 0, 1],
    #               orientation=orientation_test,
    #               velocity=[0, 0, 0, 0, 0, 0],
    #               mass=1, inertia=inertia_test)

    import math
    cs = math.cos(math.pi/4.0)
    ss = math.sin(math.pi/4.0)
    orientation_rot_x= [cs, ss, 0, 0]

    io.add_object('cyl_rot_2', [Contactor('Cyl')],
                  translation=[1, 0, L/2.0-.1],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    
    io.add_object('cyl_rot_1', [Contactor('Cyl')],
                  translation=[1, 0, L/2.0+L-0.1],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)



    io.add_object('cyl_rel_rot_1', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[0, 0, L/2.0-0.1],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    io.add_object('cyl_rel_rot_2', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[0, 0, L/2.0-0.1+L],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    
    cs = math.cos(math.pi/8.0)
    ss = math.sin(math.pi/8.0)
    orientation_rot_x= [cs, ss, 0, 0]


    io.add_object('cyl_rot_3', [Contactor('Cyl')],
                  translation=[3, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rot_4', [Contactor('Cyl')],
                  translation=[3, L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rel_rot_3', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[2, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    
    io.add_object('cyl_rel_rot_4', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[2,  L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    cs = math.cos(math.pi/4.0)
    ss = math.sin(math.pi/4.0)
    orientation_rot_x= [cs, ss, 0, 0]
    orientation_rot_z= [cs, 0, 0, ss]

    io.add_object('cyl_rot_5', [Contactor('Cyl', relative_orientation=orientation_rot_z)],
                  translation=[6, 0, 0.1],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rot_6', [Contactor('Cyl', relative_orientation=orientation_rot_z)],
                  translation=[6, 0, 0.1+2*R],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rel_rot_5', [Contactor('Cyl', relative_orientation=orientation_rot_x)],
                  translation=[4, 0, (L/2.0-.1)],
                  orientation= orientation_rot_z,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    
    io.add_object('cyl_rel_rot_6', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[4,  0, (L/2.0+L-0.1)],
                  orientation= orientation_rot_z,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)



    io.add_object('cyl_rot_7', [Contactor('Cyl')],
                  translation=[8, 0, L/2.0-.1],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    io.add_object('cyl_rot_8', [Contactor('Cyl')],
                  translation=[9, 0, L/2.0+L-0.1],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    io.add_object('cyl_rel_rot_7', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[9, 0, L/2.0-0.1],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    io.add_object('cyl_rel_rot_8', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[8, 0, L/2.0-0.1+L],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)   
    
    cs = math.cos(math.pi/8.0)
    ss = math.sin(math.pi/8.0)
    orientation_rot_x= [cs, ss, 0, 0]


    io.add_object('cyl_rot_9', [Contactor('Cyl')],
                  translation=[10, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rot_10', [Contactor('Cyl')],
                  translation=[11, L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('cyl_rel_rot_9', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[11, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)
    
    io.add_object('cyl_rel_rot_10', [Contactor('Cyl',relative_orientation=orientation_rot_x)],
                  translation=[10,  L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    
    io.add_object('compound_1', [Contactor('Cyl', relative_translation=[-0.5,0,0],  relative_orientation= orientation_rot_x),
                               Contactor('Cyl', relative_translation=[0.5,0,0], relative_orientation= orientation_rot_x)
                               ],
                  translation=[12, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('compound_rel_1', [Contactor('Cyl', relative_translation=[-0.5,0,0]),
                                   Contactor('Cyl', relative_translation=[0.5,0,0])
                               ],
                  translation=[14, 0, (L/2.0-.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('compound_2', [Contactor('Cyl', relative_translation=[-0.5,0,0],  relative_orientation= orientation_rot_x),
                               Contactor('Cyl', relative_translation=[0.5,0,0], relative_orientation= orientation_rot_x)
                               ],
                  translation=[12, L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  #orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    io.add_object('compound_rel_2', [Contactor('Cyl', relative_translation=[-0.5,0,0]),
                                   Contactor('Cyl', relative_translation=[0.5,0,0])
                               ],
                  translation=[14, L*math.cos(math.pi/4.0), (L/2.0+L-0.1)*math.cos(math.pi/4.0)],
                  orientation= orientation_rot_x,
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1, inertia=inertia_test)

    

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[9, 0, -0.7])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

test= True
if test:
    T=0.3
else:
    T=20.0


with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(with_timer=False,
           face_class=None,
           edge_class=None,
           t0=0,
           T=T,
           h=0.0005,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=None,
           output_contact_index_set=0)
