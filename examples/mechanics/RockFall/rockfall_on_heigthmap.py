#
# Example of how to use the heightmap object to interact with a static terrain.
#
import numpy as np
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner

import siconos.numerics as sn
import siconos.kernel as sk

# 1 - build the heightmap as a numpy 2D array
# The matrix consists of elevation (z) of the terrain point
# over a rectangular grid.

# Generate a numpy matrix defining a slope + sinusoid + noise
noise_amp = 1
sine_amp = 2
sine_freq = 0.3
slope_offset = 50
slope_height = 50

s = np.array([np.cos(np.linspace(0, 100, 100)*sine_freq)*sine_amp]*100)
heightmap = (np.random.uniform(0, noise_amp, size=(100, 100))
             + s + s.T
             + [np.linspace(slope_offset,
                            slope_offset+slope_height, 100)]*100)

# dimension of the terrain
heightmap_size_x = 50
heightmap_size_y = 50

# 2 - Generate various type of objects that flow on terrain under gravity


# 2.1 A sphere

sphere_size = 2
sphere_mass = 1
sphere_initial_position = [0, 18, 150] # x, y, z

# 2.2 A cube

cube_size = 10
cube_mass = 1
cube_initial_position = [0, 18, 150] # x, y, z

# 2.2 A rectangle

rectangle_l = 5
rectangle_L = 3
rectangle_H = 1

rectangle_mass = 1
rectangle_initial_position = [-5, 18, 150] # x, y, z


# 2.2 A randomconvexhull defined by its vertices.
import random
from siconos.mechanics.collision.convexhull import ConvexHull

polyhedron_size = 10.0
density = 2500
polyhedron_initial_position = [5, 18, 150] # x, y, z
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



print('polyhedron_vertices', polyhedron_vertices)
ch = ConvexHull(polyhedron_vertices)
cm = ch.centroid()
print('cm', cm)
          
# correction of vertices such that o is the centroid
# and rescale it
scale = polyhedron_size / max(np.array(polyhedron_vertices).max(axis=0)
                              - np.array(polyhedron_vertices).min(axis=0))

polyhedron_vertices = (np.array(polyhedron_vertices)[:]-cm[:])*scale


print('corrected polyhedron_vertices', polyhedron_vertices)
ch = ConvexHull(polyhedron_vertices)
cm = ch.centroid()
print('cm', cm)
          
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

    # Definition of a sphere
    io.add_primitive_shape('Ball', 'Sphere', (sphere_size,))

    # Definition of a cube
    io.add_primitive_shape('Cube', 'Box', (cube_size, cube_size, cube_size))
    
    # Definition of a rectangle
    io.add_primitive_shape('Rectangle', 'Box', (rectangle_l, rectangle_L, rectangle_L))

    # Definition of a polyhedron as a convex shape
    io.add_convex_shape('Polyhedron', polyhedron_vertices)
          

    # Definition of a non smooth law. We put the objects and heightmap
    # into different collision groups so that there are object-terrain
    # collisions but no object-object collisions.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3, e=0.0,
                                      collision_group1=0, collision_group2=1)
    
    # Add a nice big sphere for effect
    io.add_object('ball',
                  [Contactor('Ball', collision_group=1)],
                  translation=sphere_initial_position,
                  mass=sphere_mass)

    # Add a nice big cube for effect
    io.add_object('cube',
                  [Contactor('Cube', collision_group=1)],
                  translation=cube_initial_position,
                  mass=cube_mass)
    
    # Add a nice big rectangle for effect
    io.add_object('rectangle',
                  [Contactor('Rectangle', collision_group=1)],
                  translation=rectangle_initial_position,
                  mass=rectangle_mass)
    
    # Add a nice big polyhedron for effect
    obj_polyhedron=io.add_object('polyhedron', [Contactor('Polyhedron', collision_group=1)],
                                 translation=polyhedron_initial_position,
                                 orientation = polyhedron_initial_orientation,
                                 velocity=polyhedron_initial_velocity,
                                 mass=polyhedron_mass,
                                 inertia=polyhedron_inertia)
    polyhedron_id = obj_polyhedron.attrs['id']
 

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('MyTerrain', collision_group=0)],
                  translation=[0, 0, 0])

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
    
with MechanicsHdf5Runner(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.

    # We run the simulation with relatively low precision because we
    # are interested mostly just in the location of contacts with the
    # height field, but we are not evaluating the performance of
    # individual contacts here.
    io.run(with_timer=False,
           verbose=True,
           verbose_progress=True,
           t0=0,
           T=T,
           h=0.01,
           theta=0.50001,
           Newton_max_iter=1,
           solver_options=options,
           numerics_verbose=False,
           output_frequency=10)


    
with MechanicsHdf5Runner(mode='r+') as io:
    data = io._out['data'] # access to hdf5 group data
    positions=data['dynamic'][:]
    # export in numpy array
    np.savetxt('positions.txt', positions)
    #print('positions')
    #print('time', 'id', ' x, y, z, q0, q1, q2, q3')
    print('positions of the polydron')
    np.set_printoptions(precision=3)
    for e in positions:
        if e[1] == polyhedron_id:
            print('at time {0:.2f} position {1} \t\t\t orientation {2}'.format(e[0], e[2:5], e[5:10]))
    
