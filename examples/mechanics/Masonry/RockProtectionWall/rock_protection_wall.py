#!/usr/bin/env python

# import subprocess
# command = "cd /Users/vincent/siconos; git symbolic-ref --short HEAD ; git lg -1; cd - ;"
# ret = subprocess.run(command, capture_output=True, shell=True)
# print(ret.stdout.decode())


from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options
from siconos.mechanics.collision.bullet import SiconosBulletOptions
import siconos.numerics as sn
import siconos.kernel as sk
import math, random
import numpy as np
# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.


Fremond=False

T = 4.0
#T = 3e-2
h_step = 5e-4

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.
bullet_options.contactBreakingThreshold = 0.04*bullet_options.worldScale
bullet_options.perturbationIterations = 3.
bullet_options.minimumPointsPerturbationThreshold = 3.
bullet_options.extrapolationCoefficient = 0.5*h_step

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-06
options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 100



#options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_ADMM)
#options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 5000
#options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-10
#options.iparam[sn.SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = sn.SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY
run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=T
run_options['h']=h_step
run_options['theta']=0.5

run_options['constraints_activation_threshold']=1e-01
run_options['activate_with_negative_relative_velocity']=True
run_options['constraint_activation_threshold_velocity']=1e-05


run_options['Newton_max_iter'] =5
run_options['Newton_tolerance'] =1e-08

run_options['bullet_options']=bullet_options
run_options['solver_options']=options


#run_options['skip_last_update_output']=True
#run_options['skip_reset_lambdas']=True
run_options['osns_assembly_type']= sk.REDUCED_DIRECT


run_options['verbose']=True
run_options['with_timer']=True
#run_options['explode_Newton_solve']=True
#run_options['explode_computeOneStep']=True
#run_options['numerics_verbose']=True
#run_options['numerics_verbose_level']=1
#run_options['violation_verbose'] = True
run_options['output_frequency']=1


run_options['output_energy_work']=True
run_options['osi_explicit_Jacobians_of_relations']=True


#configuration = 'pyramid_wall'
configuration = 'wide_wall'
configuration = 'wide_wall_with_buttress'
configuration = 'tall_wall_with_buttress'
#configuration = 'tall_wall'
#configuration = 'one_brick'

rock_velocity = [-10,0,-35.,0.5,10.,0.1]
#rock_velocity = [-10,0,-35.,0.0,0.0,0.0]
#rock_velocity = [-0,0,-0.,0.0,0.0,0.0]
rock_tob = 0.01
rock_size = 5.0

def one_rock(io, name, cname, rock_size=0.05, density=1, trans=None, velo=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]
    rd = [1.143194140731269, 0.6247636994959747, 0.4198206059540749, 1.1116480956118107, 0.8965451614785596, 0.8819019228647785, 0.2592675582459427, 0.22899315888913663, 0.23837569282753324, 0.6585606791241505, 1.0084563758002816, 0.9096785924063177, 0.8633716705085941, 1.1215975890657788, 1.1983269522076825, 0.5443688021200721]

    # print('rd', rd)
    # input()
    def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))

    vertices = [ vert( 0,  1,   1,  1,  1),
                 vert( 2,  3,   1, -1,  1),
                 vert( 4,  5,  -1,  1,  1),
                 vert( 6,  7,  -1, -1,  1),
                 vert( 8,  9,   1,  1, -1),
                 vert(10, 11,   1, -1, -1),
                 vert(12, 13,  -1,  1, -1),
                 vert(14, 15,  -1, -1, -1) ]

    scale = rock_size / max(np.array(vertices).max(axis=0)
                            - np.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (np.array(vertices)[:] - cm[:]) * scale

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.1*rock_size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)


    return io.add_object(name,
                         [Contactor(cname)],
                         translation=trans,
                         velocity=velo,
                         mass=volume*density,
                         time_of_birth=tob,
                         inertia=inertia*density)



if Fremond:
    fn='rock_protection_wall_FremondNSL.hdf5'
else:
    fn='rock_protection_wall_NewtonNSL.hdf5'



# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner(io_filename=fn) as io:

    width, depth, height = 1, 2, 1

    margin =0.01
    io.add_primitive_shape('Box', 'Box', [width, depth, height], outsideMargin=margin)
    io.add_primitive_shape('Half_Box', 'Box', [width, depth/2.0, height], outsideMargin=margin)

    k = 0
    sep = 0.01

    def make_stack(X, Y, N, M, W):
        global k
        z = height/2.0
        while W > 0:
            for i in range(N):
                for j in range(M):
                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-M/2.0)*(depth+sep) + Y
                    io.add_object('box%03d' % k, [Contactor('Box')],
                                  translation=[x, y, z],
                                  mass=1.0)
                    k += 1
            N = N - 1 if N > 1 else N
            M = M - 1 if M > 1 else M
            W = W - 1
            z += height + sep

    def make_large_stack(X, Y, N, M, W, orientation=[1.0,0.0,0.0], density=1.0):
        """
        N is the depth of the wall
        W is th height
        """
        global k
        z = height/2.0
        volume= height*width*depth
        angle = math.acos(orientation[0]/np.linalg.norm(orientation))
        while W > 0:
            for i in range(N):
                for j in range(M):
                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-M/2.0)*(depth+sep) + Y

                    # apply rotation around z axis
                    X_r =  x*math.cos(angle) + y* math.sin(angle)
                    Y_r = -x*math.sin(angle) + y* math.cos(angle)



                    io.add_object('box%03d' % k, [Contactor('Box')],
                                  translation=[X_r, Y_r, z],
                                  orientation =[math.cos(-angle/2.0),
                                                0.0,0.0,
                                                math.sin(-angle/2.0)],
                                  mass=volume*density)
                    k += 1
            N = N - 1 if N > 1 else N
            M = M - 1 if M > 1 else M
            W = W - 1
            z += height + sep

    def make_large_stack_with_straight_side(X, Y, N, M, W, orientation=[1.0,0.0,0.0], density=1.0):
        """
        N is the depth of the wall
        M is the width
        W is the height
        """
        global k
        z = height/2.0

        angle = math.acos(orientation[0]/np.linalg.norm(orientation))
        while W > 0:
            for i in range(N):
                for j in range(M):
                    if (j == M-1 and W%2 ==0):
                        shape = 'Half_Box'
                        volume= height*width*depth/2.0
                        offset = -depth/4.0
                    elif (j == 0 and W%2 ==1):
                        shape = 'Half_Box'
                        volume= height*width*depth/2.0
                        offset = depth/4.0
                    else:
                        volume= height*width*depth
                        shape = 'Box'
                        offset=0.0



                    x = (i-N/2.0)*(width+sep) + X
                    y = (j-(M +  W%2 )/2.0)*(depth+sep) + Y +offset

                    # apply rotation around z axis
                    X_r =  x*math.cos(angle) + y* math.sin(angle)
                    Y_r = -x*math.sin(angle) + y* math.cos(angle)




                    io.add_object('box%03d' % k, [Contactor(shape)],
                                  translation=[X_r, Y_r, z],
                                  orientation =[math.cos(-angle/2.0),
                                                0.0,0.0,
                                                math.sin(-angle/2.0)],
                                  mass=volume*density)
                    k += 1
            N = N - 1 if N > 1 else N
            M = M  if M > 1 else M
            W = W - 1
            z += height + sep
    # A wall

    if configuration == 'wide_wall':
        make_large_stack(0, 0, 1, 20, 5, density=2300)
    if configuration == 'tall_wall':
        make_large_stack(0, 0, 1, 20, 10, density=2300)
    elif  configuration == 'pyramid_wall':
        make_large_stack(0, 0, 1, 10, 10, density=2300)
    elif configuration == 'wide_wall_with_buttress':
        make_large_stack(0, 0, 1, 20, 5, density=2300)
        make_large_stack_with_straight_side(10,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(5,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(0,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(-5,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(-10,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
    elif configuration == 'tall_wall_with_buttress':
        buttress_offset = -depth/2.0
        make_large_stack(0, 0, 1, 20, 10, density=2300)
        make_large_stack_with_straight_side(15+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(10+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(5+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(0+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(-5+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
        make_large_stack_with_straight_side(-10+buttress_offset,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0], density=2300)
    elif configuration == 'one_brick':
        pass
        #make_large_stack(0, 0, 1, 1, 1, density=2300)

    # # a big ball
    # io.add_primitive_shape('Ball', 'Sphere', [width])
    # io.add_object('Ball' % k, [Contactor('Ball')],
    #               translation=[0., 0., 0.],
    #               mass=1.0)


    # a big rock
    obj_projectile = one_rock(io, 'rock', 'rock_shape', rock_size = rock_size,
                             density=2300,
                             trans = [30.,0. ,35.],
                             velo=rock_velocity, tob=rock_tob)

    obj_projectile_id=obj_projectile.attrs['id']


    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (50, 50, 0.1))
    io.add_object('ground', [Contactor('Ground')], [0, 0, -0.05])
    # Definition of the slope
    angle_slope= -math.pi/4.0
    io.add_object('slope', [Contactor('Ground',
                                      relative_translation=[30.0, 0.0, -25.0* math.sin(angle_slope)],
                                      relative_orientation=[math.cos(angle_slope/2.0),
                                                            0.0,
                                                            math.sin(angle_slope/2.0),
                                                            0.0]) ],
                  [0, 0, -0.05])

    # # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[25,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    if Fremond:
        io.add_Fremond_impact_friction_nsl('contact', mu=0.6, e=0.2)
    else:
        io.add_Newton_impact_friction_nsl('contact', mu=0.6, e=0.2)




run_options['time_stepping']=None

# Load and run the simulation
with MechanicsHdf5Runner(mode='r+',io_filename=fn) as io:
    io.run(run_options)
