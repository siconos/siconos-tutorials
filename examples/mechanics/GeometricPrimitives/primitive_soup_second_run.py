#!/usr/bin/env python
from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_run import MechanicsHdf5Runner
from siconos.mechanics.collision.convexhull import ConvexHull

import siconos.numerics as sn
import siconos.kernel as sk
import numpy
# A collection of box stacks for stress-testing Siconos solver with
# chains of contacts.


options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-8

from siconos.mechanics.collision.bullet import SiconosBulletOptions
bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.0
bullet_options.contactBreakingThreshold = 1.0
bullet_options.perturbationIterations = 3
bullet_options.minimumPointsPerturbationThreshold = 3

test=True
if test:
    T=1.8
    options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 100
    options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-3
else:
    T=10.0

filename_first_run='primitive_soup_first_run'
filename_second_run='primitive_soup_second_run.hdf5'

import os
import shutil


if os.path.exists(filename_first_run +  '.hdf5'):
    shutil.copyfile(filename_first_run +  '.hdf5', filename_second_run)
else:
    print('first run file is not existing')
    exit()
    
final_time_first_run = 1.


    
# Load and run the simulation
with MechanicsHdf5Runner(mode='r+', io_filename=filename_second_run) as io:
    ######### Projectile
    rock_length = 3.5
    rock_density=2650.
    rock_tob=final_time_first_run+0.01
    rock_tod=rock_tob+0.5
    def create_vertices_rock3D_ETAG(L1,L2,L3): 
        L_b = 1.
        vertices = [[-L_b / 2, -L_b / 4, -L_b / 4], [-L_b / 2, -L_b / 4, L_b / 4], [-L_b / 2, L_b / 4, -L_b / 4], [-L_b / 2, L_b / 4, L_b / 4],
                    [L_b / 2, -L_b / 4, -L_b / 4], [L_b / 2, -L_b / 4, L_b / 4], [L_b / 2, L_b / 4, -L_b / 4], [L_b / 2, L_b / 4, L_b / 4],
                    [-L_b / 4, L_b / 2, -L_b / 4], [-L_b / 4, L_b / 2, L_b / 4], [L_b / 4, L_b / 2, -L_b / 4], [L_b / 4, L_b / 2, L_b / 4],
                    [-L_b / 4, -L_b / 2, -L_b / 4], [-L_b / 4, -L_b / 2, L_b / 4], [L_b / 4, -L_b / 2, -L_b / 4], [L_b / 4, -L_b / 2, L_b / 4],
                    [-L_b / 4, -L_b / 4, L_b / 2], [-L_b / 4, L_b / 4, L_b / 2], [L_b / 4, -L_b / 4, L_b / 2], [L_b / 4, L_b / 4, L_b / 2],
                    [-L_b / 4, -L_b / 4, -L_b / 2], [-L_b / 4, L_b / 4, -L_b / 2], [L_b / 4, -L_b / 4, -L_b / 2], [L_b / 4, L_b / 4, -L_b / 2]]
        vertices=numpy.array(vertices)
        vertices[:,0]=vertices[:,0]-numpy.min(vertices[:,0])
        vertices[:,1]=vertices[:,1]-numpy.min(vertices[:,1])
        vertices[:,2]=vertices[:,2]-numpy.min(vertices[:,2])
        vertices[:,0]=L1*vertices[:,0]/numpy.max(vertices[:,0])
        vertices[:,1]=L2*vertices[:,1]/numpy.max(vertices[:,1])
        vertices[:,2]=L3*vertices[:,2]/numpy.max(vertices[:,2])
        return vertices

    vertices=create_vertices_rock3D_ETAG(rock_length,rock_length,rock_length)
    ch = ConvexHull(vertices)
    cm = ch.centroid()
    # move the vertices to center the center of mass at 0.0
    vertices = numpy.array(vertices)[:] - cm[:]
    ch = ConvexHull(vertices)
    cm = ch.centroid()
    inertia, volume = ch.inertia(cm)

    name_new_object='aprojectile'
    #name_new_object='projectile'

    
    io.add_convex_shape(name_new_object+ '_shape', vertices,outsideMargin=0.0)
    io.add_object(name_new_object+'_DS', [Contactor(name_new_object+ '_shape')],translation=[0, 0, 8],velocity=[0., 0, -10., 0, 0, 0],
                  mass=rock_density*volume,inertia=inertia,time_of_birth=rock_tob, time_of_death=rock_tod)


    io.add_Newton_impact_friction_nsl('contact', mu=0.5)


    
    io.run(t0=0,
           T=T,
           h=0.001,
           theta=0.5,
           Newton_max_iter=2,
           solver_options=options,
           bullet_options=bullet_options,
           output_frequency=1)
