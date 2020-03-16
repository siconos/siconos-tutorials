#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as sn
import siconos.kernel as sk

import read_tess

import sys
import numpy

if (len(sys.argv) < 2):
    dist = 'uniform'
    mu = 0.1
else:
    dist = sys.argv[1]
    mu = sys.argv[2]

if not dist in ['uniform', 'double', 'exp']:
    print("dist = [uniform | double | exp]")
    sys.exit(1)
if float(mu) < 0.1 or float(mu) > 2.0:
    print("mu = [0.1 .. 2.0]")
    sys.exit(1)


fn = 'tesselation-{0}-mu-{1}.hdf5'.format(dist,mu)


test=True

if test:
    hstep=1e-3
    T = hstep*1000
    filename='n5-id1.tess'
    fn = 'tesselation.hdf5'
else:
    hstep=1e-3
    T = hstep*10000
    filename='n100-id1.tess'
    fn = 'tesselation.hdf5'

tesselation = read_tess.read_tesselation(filename)

#print(tesselation)
#input()
from siconos.mechanics.collision.tools import Contactor
from siconos.mechanics.collision.convexhull import ConvexHull

with MechanicsHdf5Runner(mode='w', io_filename=fn) as io:

    all_faces = tesselation['face']
    all_edges = tesselation['edge']
    all_vertices = tesselation['vertex']
    #print(all_faces[:])
    #input()
    for p in tesselation['polyhedron']:
        #print('p=', p)
        #create set of vertices
        faces= p[1]
        vertices=[]
        for f in faces:
            #print('f', f)
            [face]  = [x for x in all_faces if x[0] == abs(f)]
            #print(face)
            vertice = face[1]
            vertices.extend(vertice)
        vertices=list(set(vertices))
        #print('ver', vertices)
        if (len(vertices) <=  1):
            print('the number of vertices must be mroe than 0')
        vertices_coordinates=[]
        for v in vertices:
            #print(v)
            [vertex]  = [x for x in all_vertices if x[0] == abs(v)]
            #print(vertex)
            vertices_coordinates.append([vertex[1],vertex[2],vertex[3]])
        #print(vertices_coordinates)

        # Definition of a polyhedron as a convex shape
        cname = 'polyhedron' + str(p)
        ch = ConvexHull(vertices_coordinates)
        cm = ch.centroid()
        vertices_coordinates = (numpy.array(vertices_coordinates)[:] - cm[:])
        io.add_convex_shape(cname, vertices_coordinates, insideMargin=0.0)

        ch = ConvexHull(vertices_coordinates)
        inertia,volume=ch.inertia(ch.centroid())

        density=2300
        #print('geometric inertia:', inertia)
        #print('volume:', volume)
        # print('mass:', volume*density)
        # print('inertia:', inertia*density)


        name = 'polyhedron_bdy' + str(p)

        io.add_object(name,
                      [Contactor(cname)],
                      translation=cm*1.0,
                      orientation=[1.0, 0.0, 0., 0.],
                      velocity=[0.,0.,0.,0.,0.,0.],
                      mass=volume*density,
                      inertia=inertia*density)
        #input()

    io.add_primitive_shape('Ground', 'Box', (4.0, 4.0, 0.2),
                           insideMargin=0.0, outsideMargin=0.0)
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0.5, 0.5, -0.5-0.5])

    io.add_Newton_impact_friction_nsl('contact', mu=1.0, e=0.0)



# Create solver options
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
with MechanicsHdf5Runner(mode='r+', io_filename=fn) as io:
    io.run(t0=0,
           T=T,
           h=hstep,
           multipoints_iterations=True,
           theta=1.0,
           Newton_max_iter=1,
           solver_options=options,
           output_frequency=1)
