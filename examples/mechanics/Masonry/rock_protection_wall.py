#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics
import math, random, numpy
# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.

def one_rock(io, name, cname, rock_size=0.05, density=1, trans=None, velo=None, tob=None):
    # Definition of an irregular polyhedron as a convex shape

    rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]

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

    scale = rock_size / max(numpy.array(vertices).max(axis=0)
                            - numpy.array(vertices).min(axis=0))

    ch = ConvexHull(vertices)
    cm = ch.centroid()

    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm[:]) * scale

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


    io.add_object(name,
                 [Contactor(cname)],
                 translation=trans,
                 velocity=velo,
                 mass=volume*density,
                 time_of_birth=tob,
                 inertia=inertia*density)




#configuration = 'pyramid_wall'
#configuration = 'wide_wall'
configuration = 'wide_wall_with_butress'

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    width, depth, height = 1, 2, 1
    io.add_primitive_shape('Box', 'Box', [width, depth, height])
    io.add_primitive_shape('Half_Box', 'Box', [width, depth/2.0, height])

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
        angle = math.acos(orientation[0]/numpy.linalg.norm(orientation))
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

        angle = math.acos(orientation[0]/numpy.linalg.norm(orientation))
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

    if configuration == 'wide_wall' or configuration =='wide_wall_with_butress':
        make_large_stack(0, 0, 1, 20, 5)
    elif  configuration == 'pyramid_wall':
        make_large_stack(0, 0, 1, 10, 10)
        
    if configuration == 'wide_wall_with_butress':
        make_large_stack_with_straight_side(10,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0])
        make_large_stack_with_straight_side(5,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0])
        make_large_stack_with_straight_side(0,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0])
        make_large_stack_with_straight_side(-5,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0])
        make_large_stack_with_straight_side(-10,-1.5*depth, 1, 4, 5, orientation= [0.0,1.0, 0.0])
        
    # # a big ball
    # io.add_primitive_shape('Ball', 'Sphere', [width])
    # io.add_object('Ball' % k, [Contactor('Ball')],
    #               translation=[0., 0., 0.],
    #               mass=1.0)
    

    # a big rock
    one_rock(io, 'rock', 'rock_shape', rock_size = 2.5* width,
             density=2300,
             trans = [30.,0. ,35.],
             velo=[-5,0,-15.,0.5,0.1,0.1], tob=0.01 )

    
    
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
    io.add_Newton_impact_friction_nsl('contact', mu=0.6, e=0.0)

T = 6.0
#T = 1e-2
h_step = 5e-3
    
# Load and run the simulation
with MechanicsHdf5Runner(mode='r+') as io:
    io.run(t0=0,
           T=T,
           h=h_step,
           theta=0.5,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-04,
           output_frequency=1)
