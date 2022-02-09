#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
# using the Siconos proposed mechanics API
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options

import siconos.numerics as sn
import siconos.kernel as sk

import math
# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a sphere
    io.add_primitive_shape('Sphere', 'Sphere', (2,),
                           insideMargin=0.2, outsideMargin=0.0)

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (10, 10, 0.1),
                           insideMargin=0.05, outsideMargin=0.0)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.1, e=0.0)

    # The sphere object made with an unique Contactor : the sphere shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    sphere=io.add_object('sphere', [Contactor('Sphere')],
                  translation=[0, 0, 2.5],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1)
    
    io.add_object('sphere2', [Contactor('Sphere')],
                  translation=[-2.5, 0., 5.5],
                  velocity=[0, 0, 0, 0, 0, 0],
                  mass=1)
    
    # io.add_object('sphere3', [Contactor('Sphere')],
    #               translation=[1.5, 0., 7.5],
    #               velocity=[0, 0, 0, 0, 0, 0],
    #               time_of_birth=0.02,
    #               time_of_death=1.5,
    #               mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[2, 0, 0.3])

    angle = math.pi/2.0
    static_body=io.add_object('ground-tob', [Contactor('Ground')],
                              translation=[2.5, 0, -0.3],
                              orientation=[math.cos(angle/2.0), 0, math.sin(angle/2.0), 0.],
                              time_of_birth=0.001)

    static_body_id = static_body.attrs['id']
    print("static_body_id", static_body_id)


    
    angle = -math.pi/3.0
    ground_tob2=io.add_object('ground-tob-2', [Contactor('Ground')],
                              translation=[3.0, 0, -0.3],
                              orientation=[math.cos(angle/2.0), 0, math.sin(angle/2.0), 0.],
                              time_of_birth=0.55)
    

angle = math.pi/4.0
    
def apply_gravity(body):
    g = 9.81
    weight = [body.scalarMass() * g * math.sin(angle), 0.,
              - body.scalarMass() * g * math.cos(angle)]
    body.setFExtPtr(weight)  # scalMass() dans quel bibli ?

    


import numpy
from siconos.mechanics.collision.bullet import *
class death_hook():
    def __init__(self):
        pass

    def initialize(self, io):
        self._io= io
        pass

    def call(self, step):
        print('call death hook at step', step)

        # print(self._io._input.keys())
        # obj=self._io._input['ground-tob']
        # print(obj.keys())
        # print(obj.items())
        # #input_ctrs= [ctr for _n_, ctr in obj.items()])
        
        # data = obj['Ground-0']
        # print(data.attrs['instance_name'])
        # input_ctrs= [ctr for _n_, ctr in obj.items()]
        
        # input()
        
        # list interactions
        # interactions  = self._io._nsds.InteractionsVector()
        # for i in interactions:
        #     print("interaction number", i.number())
        #     #i.display()

        #list contact points
        contact_points = self._io._io.contactPoints(self._io._nsds, 1)
        if contact_points is not None:
            #print('contact_points')
            #print(contact_points )
            for cp in contact_points:
                #print("cp", cp)
                inter_id = int(cp[22])
                #print("inter_id", inter_id)
                ds1_id=int(cp[23])
                ds2_id=int(cp[24])
                #print("ds id :", ds1_id, ds2_id)
                if ds1_id == ds2_id:
                    print("interaction with a static object")
                    inter = self._io._nsds.interaction(inter_id)
                    
                    contact_r = cast_BulletR(inter.relation())
                    print("contact_r.bodyShapeRecordA", contact_r.bodyShapeRecordA)
                    print("contact_r.bodyShapeRecordB", contact_r.bodyShapeRecordB)
                    
                    print("contact_r.bodyShapeRecordB.staticBody")
                    print("contact_r.bodyShapeRecordB.staticBody", contact_r.bodyShapeRecordB.staticBody)
                    print("contact_r.bodyShapeRecordB.staticBody.number", contact_r.bodyShapeRecordB.staticBody.number)
                    print("contact_r.bodyShapeRecordB.ds", contact_r.bodyShapeRecordB.ds)
                    
                    print("contact_r.bodyShapeRecordA.staticBody")
                    print("contact_r.bodyShapeRecordA.staticBody", contact_r.bodyShapeRecordA.staticBody)
                    print("contact_r.bodyShapeRecordA.staticBody.ds.number()", contact_r.bodyShapeRecordA.ds.number())
                    
                    if contact_r.bodyShapeRecordB.staticBody.number == static_body_id:
                        self._io._interman.removeStaticBody(contact_r.bodyShapeRecordB.staticBody)
                        # remove the body from that list of static object completely
                        for s in self._io._static:
                            #print(self._io._static[s])
                            #print(self._io._static[s]['number'])
                            if self._io._static[s]['number'] == static_body_id:
                                s_remove=s
                        self._io._static.pop(s_remove)
        else :
            #print('no contact points')
            pass
        
        # second way (faster) :  direct access to nsds positions
        # positions  = self._io._io.positions(self._io._nsds)
        # if positions is not None:
        #     z = positions[:,3]
        #     # We search for the ds index that are below a given criteria
        #     ds_idx = numpy.nonzero(z < -2)[0]
        #     for i in ds_idx :
        #         n_ds = int(positions[i,0])
        #         ds = self._io._nsds.dynamicalSystem(n_ds)
        #         self._io._interman.removeBody(ds)
        #         self._io._nsds.removeDynamicalSystem(ds)
        #         print('remove ds number ', ds.number(), ' with height = ', z[i])



dh = death_hook()

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 1.0
bullet_options.perturbationIterations = 3.
bullet_options.minimumPointsPerturbationThreshold = 3.

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=3e-1
run_options['h']=1e-3
run_options['theta']=0.5

run_options['bullet_options']=bullet_options
run_options['solver_options']=options
run_options['constraint_activation_threshold']=1e-05

#run_options['start_run_iteration_hook']=sh
run_options['end_run_iteration_hook'] = dh

run_options['Newton_options']=sk.SICONOS_TS_LINEAR

# run_options['skip_last_update_output']=True
# run_options['skip_reset_lambdas']=True
# run_options['osns_assembly_type']= sk.REDUCED_DIRECT

# run_options['osns_assembly_type']= sk.GLOBAL_REDUCED
# run_options['osi']= sk.MoreauJeanGOSI
# #run_options['skip_last_update_input']=True

# run_options['verbose']=True
# run_options['with_timer']=True
# run_options['explode_Newton_solve']=True
# run_options['explode_computeOneStep']=True

#run_options['output_frequency']=None
run_options['output_contact_index_set']=0
#run_options['time_stepping']=None

run_options['set_external_forces']=None


with MechanicsHdf5Runner(mode='r+', set_external_forces=apply_gravity) as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    io.run(run_options)
