#!/usr/bin/env python

#
# Watch Escapment mechanism with opencascade contactors
#
# In order to compile the needed plugins before the execution of this
# example, use the siconos command:
#
# siconos <this file>
#

from pathlib import Path
from siconos.mechanics.collision.tools import Contactor, Shape, Volume, Material
from siconos.io.mechanics_run import MechanicsHdf5Runner
from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder
from OCC.gp import gp_Pnt, gp_Ax2, gp_Dir
import siconos.numerics as sn
import siconos.kernel as sk
import siconos.io.mechanics_run
siconos.io.mechanics_run.set_backend('occ')
import numpy as np
import math

scale = 1e-9
density = 7750
steel = Material(density=density)
steel_mm = Material(density=density*scale)

restart = False
if not restart:
    with MechanicsHdf5Runner() as io:
        step_dir = Path(__file__).resolve().parent
        step_dir = Path(step_dir, 'step_files')
        step_file = Path(step_dir, 'BalanceWheel.step').as_posix()
        io.add_shape_data_from_file('balance_wheel_shp', step_file)

        step_file = Path(step_dir, 'BalanceWheelBearings.step').as_posix()
        io.add_shape_data_from_file('balance_wheel_bearings_shp', step_file)

        step_file = Path(step_dir, 'BalanceWheelContactor.step').as_posix()
        io.add_shape_data_from_file('balance_wheel_contactor_shp', step_file)

        step_file = Path(step_dir, 'PinLever.step').as_posix()
        io.add_shape_data_from_file('pin_lever_shp', step_file)

        step_file = Path(step_dir, 'PinLeverBearings.step').as_posix()
        io.add_shape_data_from_file('pin_lever_bearings_shp', step_file)

        # contactors between Pin Lever and Escape Wheel
        step_file = Path(step_dir, 'PinLeverContactorA.step').as_posix()
        io.add_shape_data_from_file('pin_lever_contactorA_shp', step_file)

        step_file = Path(step_dir, 'PinLeverContactorB.step').as_posix()
        io.add_shape_data_from_file('pin_lever_contactorB_shp', step_file)

        # contactors between Pin Lever and back Cover (stops for Pin Lever)
        use_cylinder_as_stops = True
        if (use_cylinder_as_stops):
            step_file = Path(step_dir, 'PinLeverStopsA.step').as_posix()
            io.add_shape_data_from_file('pin_lever_stopsA_shp', step_file)
            step_file = Path(step_dir, 'PinLeverStopsB.step').as_posix()
            io.add_shape_data_from_file('pin_lever_stopsB_shp', step_file)
        else:
            # back cover step causes issues in vview
            # io.add_shape_data_from_file('back_cover_shp',
            #                             'BackCover.step')
            step_file = Path(step_dir,
                             'BackCover_Contactor_face21.step').as_posix()
            io.add_shape_data_from_file('back_cover_contactorA_shp', step_file)
            step_file = Path(step_dir,
                             'BackCover_Contactor_face213.step').as_posix()
            io.add_shape_data_from_file('back_cover_contactorB_shp', step_file)

        step_file = Path(step_dir, 'EscapeWheel.step').as_posix()
        io.add_shape_data_from_file('escape_wheel_shp', step_file)

        step_file = Path(step_dir, 'EscapeWheelBearings.step').as_posix()
        io.add_shape_data_from_file('escape_wheel_bearing_shp', step_file)

        ###############
        # Balance Wheel
        ###############

        # to compute the center of mass and perform the right translation
        # io.add_object('balance_wheel_body_com',
        #               [Volume(shape_name='balance_wheel_shp',
        #                       instance_name='balance_wheel_shp_bdy')],
        #               translation=[0., 0., 0.],
        #               velocity=[0., 0., 0., 0., 0., 0.],
        #               mass = 100.0)
        balance_wheel_com = np.array([ 10.36290606,  -3.90431343,  -0.21335326])


        inertia = [.019661740055950054, 0.028171281786527207, 0.021165152882136887]
        init_orientation = [0.97366047, 0., 0., -0.22800273] #

        angle_init = math.asin(-0.22800273)
        print('angle_init =', angle_init)
        angle_init = math.pi/4.0
        init_orientation =  [math.sin(angle_init), 0., 0., math.cos(angle_init)]

        io.add_object('balance_wheel_body',
                      [Volume(shape_name='balance_wheel_shp',
                              instance_name='balance_wheel_shp_bdy',
                              relative_translation=-1.0*balance_wheel_com,
                              parameters=steel_mm),
                       Contactor(
                           instance_name='balance_wheel_contactor_shp_f0',
                           shape_name='balance_wheel_contactor_shp',
                           contact_type='Face',
                           contact_index=0,
                           relative_translation=-1.0*balance_wheel_com),
                       Contactor(
                           instance_name='balance_wheel_contactor_shp_f1',
                           shape_name='balance_wheel_contactor_shp',
                           contact_type='Face',
                           contact_index=1,
                           relative_translation=-1.0*balance_wheel_com)],
                      translation=balance_wheel_com,
                      orientation = init_orientation,
                      velocity=[0., 0., 0., 0., 0., 0.])

        io.add_joint('balance_wheel_joint',  'balance_wheel_body',
                     points=[balance_wheel_com],
                     axes=[[0., 0., 1.]],
                     joint_class='PivotJointR',
                     absolute=True)

        io.add_external_function('f1', 'balance_wheel_body', 'setComputeMIntFunction',
                                 'Plugin', 'internalMomentsBalanceWheel')
        io.add_external_function('f1_jacq', 'balance_wheel_body', 'setComputeJacobianMIntqFunction',
                                 'Plugin', 'internalMomentsBalanceWheel_Jacq')

        # a static object (mass=0)
        io.add_object('balance_wheel_bearings_static',
                      [Shape(
                          instance_name='balance_wheel_bearings_shp',
                          shape_name='balance_wheel_bearings_shp',
                          relative_translation=[0, 0, 0])],
                      translation=[0, 0, 0])

        ###############
        # Pin Lever
        ###############

        # to compute the center of mass
        # pin_lever_shapes= [Volume(shape_name='pin_lever_shp',
        #                      instance_name='pin_lever_shp_bdy',
        #                      relative_translation=[0., 0., 0.],
        #                      relative_orientation=[(0, 1, 0), 0.])]
        # io.add_object('pin_lever_body_com',
        #               pin_lever_shapes,
        #               translation=[0., 0., 0.],
        #               velocity=[0., 0., 0., 0., 0., 0.],
        #               mass = 1.0)

        pin_lever_com = np.array([7.70341513774, -7.40214527376, -1.23418992515])

        # contactor with Escape Wheel
        pin_lever_shapes= [Volume(shape_name='pin_lever_shp',
                                  instance_name='pin_lever_shp_bdy',
                                  relative_translation=-1.0*pin_lever_com,
                                  relative_orientation=[(0, 1, 0), 0.],
                                  parameters=steel_mm)]

        pin_lever_shapes.extend([Contactor(instance_name='pin_lever_contactorA_shp_f1',
                                           shape_name='pin_lever_contactorA_shp',
                                           contact_type='Face',
                                           contact_index=1,
                                           relative_translation=-1.0*pin_lever_com),
                                 Contactor(instance_name='pin_lever_contactorB_shp_f0',
                                           shape_name='pin_lever_contactorB_shp',
                                           contact_type='Face',
                                           contact_index=0,
                                           relative_translation=-1.0*pin_lever_com)])
        # contactor with Balance Wheel 
        pin_lever_shapes.extend([Contactor(instance_name='pin_lever_shp_f45',
                                           shape_name='pin_lever_shp',
                                           contact_type='Face',
                                           contact_index=45,
                                           relative_translation=-1.0*pin_lever_com),
                                 Contactor(instance_name='pin_lever_shp_f43',
                                           shape_name='pin_lever_shp',
                                           contact_type='Face',
                                           contact_index=43,
                                           relative_translation=-1.0*pin_lever_com)])
        # contactor with Stops or Back cover
        pin_lever_shapes.extend([Contactor(instance_name='pin_lever_shp_f39',
                                           shape_name='pin_lever_shp',
                                           contact_type='Face',
                                           contact_index=39,
                                           relative_translation=-1.0*pin_lever_com),
                                 Contactor(instance_name='pin_lever_shp_f12',
                                           shape_name='pin_lever_shp',
                                           contact_type='Face',
                                           contact_index=12,
                                           relative_translation=-1.0*pin_lever_com)])

        io.add_object('pin_lever_body',
                      pin_lever_shapes,
                      translation=pin_lever_com,
                      velocity=[0., 0., 0., 0., 0., 0.])

        pin_lever_bearings_com=[7.08439152004, -7.80542761001, -1.97034359443]
        io.add_joint('pin_lever_joint',  'pin_lever_body',
                     points=[pin_lever_bearings_com],
                     axes=[[0., 0., 1.]],
                      joint_class='PivotJointR',
                     absolute=True)

        io.add_object('pin_lever_bearings_static',
                      [Shape(
                          instance_name='pin_lever_bearings_shp',
                          shape_name='pin_lever_bearings_shp',
                          relative_translation=[0, 0, 0])],
                      translation=[0, 0, 0])
        if use_cylinder_as_stops:
            io.add_object('pin_lever_stops_static',
                          [Contactor(
                              instance_name='pin_lever_stopsA_shp_f0',
                              shape_name='pin_lever_stopsA_shp',
                              contact_type='Face',
                              contact_index=0),
                           Contactor(
                               instance_name='pin_lever_stopsB_shp_f0',
                               shape_name='pin_lever_stopsB_shp',
                               contact_type='Face',
                               contact_index=0)
                          ],
                          translation=[0, 0, 0])
        else:
            io.add_object('back_cover_static',
                          [Contactor(
                              instance_name='back_cover_contactor_shp_f212',
                              shape_name='back_cover_contactorA_shp',
                              contact_type='Face',
                              contact_index=0),
                           Contactor(
                               instance_name='back_cover_contactor_shp_f213',
                               shape_name='back_cover_contactorB_shp',
                               contact_type='Face',
                               contact_index=0)
                          ],
                          translation=[0, 0, 0],
            )



        #########################
        #       Escape Wheel    #
        #########################

        # escape_wheel_shapes= [Volume(shape_name='escape_wheel_shp',
        #                      instance_name='escape_wheel_shp_bdy',
        #                      relative_translation=[0.,0.,0.] ,
        #                      relative_orientation=[(0, 1, 0), 0.])]

        # io.add_object('escape_wheel_body_com',
        #               escape_wheel_shapes,
        #               translation=[0.,0.,0.],
        #               velocity=[0., 0., 0., 0., 0., 0.],
        #               mass = 1.0)

        escape_wheel_com = np.array([4.69701693231, -10.6460029586, -0.836711176511])

        list_contactor_faces=[102, 109, 116, 123, 130, 137, 144, 151, 158, 165, 172, 74, 81, 88, 95]
        #list_contactor_faces=[]
        list_contactor=[]
        for idx in list_contactor_faces:
            list_contactor.append(Contactor(
                instance_name='escape_wheel_shp_f'+str(idx),
                shape_name='escape_wheel_shp',
                contact_type='Face',
                contact_index=idx,
                relative_translation=-1.0*escape_wheel_com))

        print(list_contactor)


        escape_wheel_shapes= [Volume(shape_name='escape_wheel_shp',
                                     instance_name='escape_wheel_shp_bdy',
                                     relative_translation=-1.0*escape_wheel_com ,
                                     relative_orientation=[(0, 1, 0), 0.],
                                     parameters=steel_mm)]

        escape_wheel_shapes.extend(list_contactor)


        io.add_object('escape_wheel_body',
                      escape_wheel_shapes,
                      translation=escape_wheel_com ,
                      velocity=[0., 0., 0., 0., 0., 0.])

        io.add_external_function('M1', 'escape_wheel_body', 'setComputeMExtFunction',
                                 'Plugin', 'externalMomentEscapeWheel')

        io.add_joint('escape_wheel_joint',  'escape_wheel_body',
                     points=[escape_wheel_com],
                     axes=[[0., 0., 1.]],
                      joint_class='PivotJointR',
                     absolute=True)

        io.add_object('escape_wheel_bearing_static',
                      [Shape(
                          instance_name='escape_wheel_bearing_shp',
                          shape_name='escape_wheel_bearing_shp',
                          relative_translation=[0, 0, 0])],
                      translation=[0, 0, 0])


        offset=0.01
        io.add_interaction('contact_balance_wheel_pin_lever_10',
                           body1_name='balance_wheel_body', contactor1_name='balance_wheel_contactor_shp_f0',
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_shp_f45',
                           distance_calculator='cadmbtb',
                           offset1=offset)

        io.add_interaction('contact_balance_wheel_pin_lever_11',
                           body1_name='balance_wheel_body', contactor1_name='balance_wheel_contactor_shp_f0',
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_shp_f43',
                           distance_calculator='cadmbtb',
                           offset1=offset)

        io.add_interaction('contact_balance_wheel_pin_lever_20',
                           body1_name='balance_wheel_body', contactor1_name='balance_wheel_contactor_shp_f1',
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_shp_f45',
                           distance_calculator='cadmbtb',
                           offset1=offset)
        io.add_interaction('contact_balance_wheel_pin_lever_21',
                           body1_name='balance_wheel_body', contactor1_name='balance_wheel_contactor_shp_f1',
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_shp_f43',
                           distance_calculator='cadmbtb',
                           offset1=offset)

        counter=0
        for c in list_contactor:

            io.add_interaction('contact_escape_wheel_pin_lever_'+str(counter),
                           body1_name='escape_wheel_body', contactor1_name=c.instance_name,
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_contactorA_shp_f1',
                           distance_calculator='cadmbtb',
                           offset1=offset)
            counter=counter+1
            io.add_interaction('contact_escape_wheel_pin_lever_'+str(counter),
                           body1_name='escape_wheel_body', contactor1_name=c.instance_name,
                           body2_name='pin_lever_body',
                           contactor2_name='pin_lever_contactorB_shp_f0',
                           distance_calculator='cadmbtb',
                           offset1=offset)
            counter=counter+1



        if use_cylinder_as_stops:
            offset=0.01
            io.add_interaction('contact_pin_lever_stops_0',
                               body1_name='pin_lever_body', contactor1_name='pin_lever_shp_f39',
                               body2_name='pin_lever_stops_static', contactor2_name='pin_lever_stopsA_shp_f0',
                               distance_calculator='cadmbtb',
                               offset1=offset)

            io.add_interaction('contact_pin_lever_stops_1',
                               body1_name='pin_lever_body', contactor1_name='pin_lever_shp_f12',
                               body2_name='pin_lever_stops_static', contactor2_name='pin_lever_stopsB_shp_f0',
                               distance_calculator='cadmbtb',
                               offset1=offset)
        else:
            offset=0.05
            io.add_interaction('contact_back_cover_pin_lever_0',
                               body1_name='back_cover_static', contactor1_name='back_cover_contactor_shp_f212',
                               body2_name='pin_lever_body',
                               contactor2_name='pin_lever_shp_f12',
                               distance_calculator='cadmbtb',
                               offset1=offset)
            io.add_interaction('contact_back_cover_pin_lever_1',
                               body1_name='back_cover_static', contactor1_name='back_cover_contactor_shp_f213',
                               body2_name='pin_lever_body',
                               contactor2_name='pin_lever_shp_f39',
                               distance_calculator='cadmbtb',
                               offset1=offset)


        io.add_Newton_impact_friction_nsl('contact', mu=0.01, e=0.0)

h_step = 1e-4
n_step = 300000

options = sk.solver_options_create(sn.SICONOS_GENERIC_MECHANICAL_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 10000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-4
sk.solver_options_update_internal(options, 1, sn.SICONOS_FRICTION_3D_ONECONTACT_NSN)


with MechanicsHdf5Runner(mode='r+') as io:

    io.run(with_timer=True,
           t0=0,
           T=n_step*h_step,
           h=h_step,
           solver_options=options,
           Newton_max_iter=5,
           contact_index_set=0,
           output_frequency=10)
