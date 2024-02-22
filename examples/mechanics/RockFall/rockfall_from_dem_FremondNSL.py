#!/usr/bin/env python

#
# Example of how to use the heightmap object to interact with a static terrain.
#
import numpy as np

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner, MechanicsHdf5Runner_run_options
from siconos.mechanics.collision.bullet import SiconosBulletOptions

import siconos.numerics as sn
import siconos.kernel as sk



# 1 - build the heightmap as a numpy 2D array
# The matrix consists of elevation (z) of the terrain point
# over a rectangular grid.

import terrain
Terrain = terrain.terrain('./data/dem.asc')
Terrain.change_min(1000)
#Terrain.rotation(-60)
print(Terrain, Terrain.dimension)
#Terrain.plot2d()
heightmap = Terrain.matrix
heightmap_size_x = Terrain.dimension[0]
heightmap_size_y = Terrain.dimension[1]

# 2 - Generate a random convexhull defined by its vertices
# that flow on terrain under gravity given by the file ./data/d1/asc

import random
from siconos.mechanics.collision.convexhull import ConvexHull
polyhedron_size = 10.0
density = 2500

polyhedron_initial_position=np.array(Terrain.initial_position('./data/d1.asc'))
polyhedron_initial_position[2] = polyhedron_initial_position[2] +  polyhedron_size
#print(polyhedron_initial_position)
import math
angle =math.pi/2.0
polyhedron_initial_orientation =   [math.cos(angle/2.0), 0.0, math.sin(angle/2.0), 0.0] # quaternion (rotation aroung the y-axis of angle)
polyhedron_initial_velocity = [0., 0., 0., 0., 0., 0.]

rd = [math.pi/2 * random.gauss(0.5,0.2) for _ in range(16)]


def vert(id1, id2, a, b, c):
        return (a*math.cos(rd[id1])*math.cos(rd[id2]),
                b*math.sin(rd[id1])*math.cos(rd[id2]),
                c*math.sin(rd[id2]))
polyhedron_vertices = [ vert( 0,  1,   1,  1,  1),
                        vert( 2,  3,   1, -1,  1),
                        vert( 4,  5,  -1,  1,  1),
                        vert( 6,  7,  -1, -1,  1),
                        vert( 8,  9,   1,  1, -1),
                        vert(10, 11,   1, -1, -1),
                        vert(12, 13,  -1,  1, -1),
                        vert(14, 15,  -1, -1, -1) ]

ch = ConvexHull(polyhedron_vertices)
cm = ch.centroid()
print('orginal centroid', cm)
          
# correction of vertices such that o is the centroid
# and rescale it
scale = polyhedron_size / max(np.array(polyhedron_vertices).max(axis=0)
                              - np.array(polyhedron_vertices).min(axis=0))

polyhedron_vertices = (np.array(polyhedron_vertices)[:]-cm[:])*scale


print('polyhedron_vertices', polyhedron_vertices)
ch = ConvexHull(polyhedron_vertices)
cm = ch.centroid()
print('centroid:', cm)
          
# computation of inertia and volume
inertia,volume=ch.inertia(cm)
print ('inertia,volume',inertia,volume)
polyhedron_mass = volume*density
polyhedron_inertia = inertia*density
print ('polyhedron inertia, polyhedron mass',polyhedron_inertia, polyhedron_mass)

# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:


    # The heightmap has tangential extents 50 by 50, and its height is
    # not scaled; actual height must be reflected by values in
    # heightmap data.
    io.add_height_map('MyTerrain', heightmap, (heightmap_size_x, heightmap_size_y))

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape('Polyhedron', polyhedron_vertices)
          
    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.add_Fremond_impact_friction_nsl('contact', mu=0.3, e=0.0,
                                      collision_group1=0, collision_group2=1)
    
    # Add a  polyhedron 
    obj_polyhedron=io.add_object('polyhedron', [Contactor('Polyhedron', collision_group=1)],
                                 translation=polyhedron_initial_position,
                                 orientation = polyhedron_initial_orientation,
                                 velocity=polyhedron_initial_velocity,
                                 mass=polyhedron_mass,
                                 inertia=polyhedron_inertia)
    polyhedron_id = obj_polyhedron.attrs['id']
    #obj_polyhedron.attrs['color'] = [0., 1.0, 0.]
 

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    obj_ground=io.add_object('ground', [Contactor('MyTerrain', collision_group=0)],
                             translation=[0, 0, 0])
    obj_ground.attrs['color'] = [0.0, 0.0, 1.0]

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-5


test=True
if test:
    T=10.0
else:
    T=30.0

time_step = 1e-2
bullet_options = SiconosBulletOptions()
# bullet_options.worldScale = .1
bullet_options.contactBreakingThreshold = 0.4
# bullet_options.perturbationIterations = 3.
# bullet_options.minimumPointsPerturbationThreshold = 3.


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-5


run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=T
run_options['h']=time_step

run_options['bullet_options']=bullet_options
run_options['solver_options']=options


#run_options['skip_last_update_output']=True
#run_options['skip_reset_lambdas']=True
#run_options['osns_assembly_type']= sk.REDUCED_DIRECT


run_options['verbose']=True
run_options['with_timer']=True
#run_options['explode_Newton_solve']=True
#run_options['explode_computeOneStep']=True

#run_options['violation_verbose'] = True
run_options['output_frequency']=1

run_options['time_stepping']=None
    
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.

    # We run the simulation with relatively low precision because we
    # are interested mostly just in the location of contacts with the
    # height field, but we are not evaluating the performance of
    # individual contacts here.
    io.run(run_options)


    
with MechanicsHdf5Runner(mode='r+') as io:
    data = io._out['data'] # access to hdf5 group data
    positions=data['dynamic'][:]
    # export in numpy array
    np.savetxt('positions.txt', positions)
    print('positions', positions)
    #print('time', 'id', ' x, y, z, q0, q1, q2, q3')
    print('positions of the polydron')
    np.set_printoptions(precision=3)
    for e in positions:
        if e[1] == polyhedron_id:
            print('at time {0:.2f} position {1} \t\t\t orientation {2}'.format(e[0], e[2:5], e[5:10]))
    
