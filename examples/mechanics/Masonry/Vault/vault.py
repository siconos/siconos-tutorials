#!/usr/bin/env python

from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull
from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics
import math, random, numpy




import pickle
brick = list(pickle.load(open( 'brick.dat', "rb" )))
vertices = list(pickle.load(open( 'vertices.dat', "rb" )))

def one_brick(io, name, cname, vertices, size, density=1, trans=None, velo=None, tob=None):


    estimated_size = max(numpy.array(vertices).max(axis=0)
                         - numpy.array(vertices).min(axis=0))
    #print(estimated_size)
    # scale = size / max(numpy.array(vertices).max(axis=0)
    #                         - numpy.array(vertices).min(axis=0))
    scale=1.0
    ch = ConvexHull(vertices)
    cm_ori = ch.centroid()
    #print('cm_ori', cm_ori)
    if (numpy.linalg.norm(cm_ori) <= 1e-2):
        input()
    # correction of vertices such that 0 is the centroid
    vertices = (numpy.array(vertices)[:] - cm_ori[:]) * scale

    
    ch = ConvexHull(vertices)
    cm = ch.centroid()
    #print('cm', cm)
    # Definition of a polyhedron as a convex shape
    io.add_convex_shape(cname, vertices, insideMargin=0.001*size)

    # computation of inertia and volume
    inertia,volume=ch.inertia(ch.centroid())

    # print('geometric inertia:', inertia)
    # print('volume:', volume)
    # print('mass:', volume*density)
    # print('inertia:', inertia*density)


    # io.add_object(name,
    #              [Contactor(cname, relative_translation = cm_ori)],
    #              translation=-cm_ori,
    #              velocity=velo,
    #              mass=volume*density,
    #              time_of_birth=tob,
    #              inertia=inertia*density)
    io.add_object(name,
                  [Contactor(cname, relative_orientation = [1.0, 0.0, 0.0, 0.0],)],
                  translation=cm_ori ,
                  orientation = [1.0, 0.0, 0.0, 0.0],
                  velocity=velo,
                  mass=volume*density,
                  time_of_birth=tob,
                  inertia=inertia*density)
    

   




# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    k=0

    vertices_array ={}
   
    i=0
    for v in vertices:
        #print(v['number'])
        vertices_array[v['number']] =i
        i=i+1
    #print(vertices_array)
    #print(brick)
    for b in brick :
        v = []
        for vb in b:
            #print('vb',vb)
            #print(len(vertices))
            v.append(vertices[vertices_array[vb]]['coord'])
        #print('v',v)
        name = 'brick%03d' % k
        cname= 'brick_shp%03d' % k
        one_brick(io, name, cname  , v, 1.0, density=2300,
                  trans=[0.0,0.0,0.0],
                  velo=[0.0,0.0,0.0, 0.0, 0.0, 0.0],tob=0.0)
        k=k+1
        # if (k > 20):
        #     break
        #input()
    
    # Definition of the ground
    io.add_primitive_shape('Ground', 'Box', (10, 1, 0.1))
    io.add_object('ground', [Contactor('Ground')], [2.5, 0, -0.5])

    # # Enable to smash the wall
    # io.add_primitive_shape('Ball', 'Sphere', [1,])
    # io.add_object('WreckingBall', [Contactor('Ball')],
    #              translation=[25,0,3], velocity=[-30,0,2,0,0,0],
    #              mass=10)

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.6, e=0.0)

T = 3.0
#T = 3e-2
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
           output_frequency=1,
           with_timer=True)
