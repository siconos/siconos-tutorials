#!/usr/bin/env python

import subprocess
command = "cd /Users/vincent/siconos; git symbolic-ref --short HEAD ; git lg -1; cd - ;"
ret = subprocess.run(command, capture_output=True, shell=True)
print(ret.stdout.decode())


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


Fremond=True

T = 2.0
#T = 3e-2
h_step = 5e-3

bullet_options = SiconosBulletOptions()
bullet_options.worldScale = 1.
bullet_options.contactBreakingThreshold = 0.04*bullet_options.worldScale
bullet_options.perturbationIterations = 3.
bullet_options.minimumPointsPerturbationThreshold = 3.

options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_NSGS)
options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 1000
options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-04
#options.iparam[sn.SICONOS_FRICTION_3D_NSGS_FREEZING_CONTACT] = 100



#options = sk.solver_options_create(sn.SICONOS_FRICTION_3D_ADMM)
#options.iparam[sn.SICONOS_IPARAM_MAX_ITER] = 5000
#options.dparam[sn.SICONOS_DPARAM_TOL] = 1e-10
#options.iparam[sn.SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY] = sn.SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY
run_options=MechanicsHdf5Runner_run_options()
run_options['t0']=0
run_options['T']=T
run_options['h']=h_step
run_options['theta']=0.5

#run_options['constraints_activation_threshold']=-1e-01
run_options['activate_with_negative_relative_velocity']=True
run_options['constraint_activation_threshold_velocity']=1e-03


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

    margin =0.
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

class kinetic_iteration_hook():

    def __init__(self, obj_projectile_id):
        self._io= None
        self._kinetic_sum = []
        self._kinetic_sum_wall = []
        self._normal_work_sum = []
        self._potential_sum = []
        self._external_forces_work = []
        self._tangent_work_sum = []
        self._friction_work_sum = []
        self._normal_work_negative_sum = []
        self._tangent_work_negative_sum = []
        self._total_energy = 0.
        self._obj_projectile_id = obj_projectile_id
        #self._frequency = 1/self._period
        #print('init hook')
        #input()
        pass

    def initialize(self, io):
        self._io= io
        pass

    def call(self, step):
        #print(' end hook step', step)
        #nsds = self._io._nsds
        positions  = self._io._io.positions(self._io._nsds)
        #velocities  = self._io._io.velocities(self._io._nsds)
        nsds= self._io._nsds


        kinetic_sum = 0.
        kinetic_sum_wall = 0.
        potential_sum = 0.
        external_forces_work_sum = 0.
        gen_contact_work=0.
        kinetic_sum_old = 0.
        time= self._io._simulation.startingTime()


        work_forces=self._io._osi.computeWorkForces()
        if work_forces is not None:
            #print(work_forces)
            external_forces_work_sum = np.sum(work_forces[:,1])



        external_forces_work_sum_adhoc = 0.
        if positions is not None:
            ds_idx  = positions[:,0]
            #print(ds_idx)

            for i in ds_idx :
                n_ds = int(i)
                ds = nsds.dynamicalSystem(n_ds)
                #ds.display()

                neds= sk.cast_NewtonEulerDS(ds)
                kinetic = neds.computeKineticEnergy()

                #print('kinetic', kinetic, 'time ', self._io._simulation.startingTime())
                if n_ds == self._obj_projectile_id:
                    #print('\nkinetic obj {0:3e}'.format(kinetic))
                    pass
                else:
                    kinetic_sum_wall =  kinetic_sum_wall+kinetic

                kinetic_sum = kinetic + kinetic_sum
                m = neds.mass()


                vold = np.array(neds.velocityMemory().getSiconosVector(0))
                v= neds.velocity()

                kinetic_old = 0.5 * np.dot(vold,np.dot(m,vold)) # compute with new mass
                kinetic_sum_old = kinetic_old + kinetic_sum_old


                pos_z = neds.q()[2]
                mass=  m[0,0]
                potential =  mass * 9.81 * pos_z
                #print('potential', potential, 'time ', self._io._simulation.startingTime())
                potential_sum = potential  + potential_sum


                #print('vold', vold)

                external_forces_work =  - h_step * 0.5 *(v[2]+ vold[2])* mass * 9.81
                external_forces_work_sum_adhoc = external_forces_work_sum+ external_forces_work

                #compute residu
                R = np.dot(m, neds.velocity()- vold)

                qold = neds.qMemory().getSiconosVector(0)
                neds.computeForces(0.0, qold, vold)
                fold = neds.forces().copy()

                neds.computeForces(0.0, neds.q(), neds.twist())
                f = neds.forces()
                # print('qold', qold)
                # print('q', ds.q())
                # print('vold', vold)
                # print('v', neds.twist())
                # print('fold',fold)
                # print('f',f)
                R  = R - h_step*0.5*(f+fold)



                R  = R - neds.p(1)

                v_k_theta = 0.5 *(v+vold)


                gen_contact_work=np.dot(v_k_theta,neds.p(1))



                #print('neds.p(1)', neds.p(1))
                # print('force work', h_step*np.dot(v_k_theta,0.5*(f+fold)) )
                # print('contact work', np.dot(v_k_theta,neds.p(1)))

                # print('norm R', np.linalg.norm(R))






        self._kinetic_sum.append([time,kinetic_sum])
        self._kinetic_sum_wall.append([time,kinetic_sum_wall])
        self._potential_sum.append([time,potential_sum])
        # print('kinetic_sum {0:3e}'.format(kinetic_sum))
        # print('potential_sum {0:3e}'.format(potential_sum))


        # print('external_forces_work_sum_adhoc', external_forces_work_sum_adhoc)
        # print('external_forces_work_sum', external_forces_work_sum)
        # print('diff', external_forces_work_sum-external_forces_work_sum_adhoc)


        self._external_forces_work.append(external_forces_work_sum)
        external_forces_work_cum = np.cumsum(self._external_forces_work)

        # print('external_forces_work_cum {0:3e}'.format(external_forces_work_cum[-1]))
        # print('diff {0:3e}'.format(external_forces_work_cum[-1]- potential_sum))

        external_forces_work_cum = np.cumsum(self._external_forces_work)

        cf_work = self._io._io.contactContactWork(self._io._nsds,
                                                  1)
        cf = self._io._io.contactPoints(self._io._nsds,
                                                  1)
        
        if cf_work is not None:
            #print('cf_work', cf_work)
            #print('cf', cf[:,22])
            interactions_idx= cf_work[:,0]
            
            contact_force=[]
            work=[]
            work_new=[]
            work_gen=0.

            P = np.zeros(6)
            #idx = 0
            for i in interactions_idx:

                i_interaction = int(i)
                inter = nsds.interaction(i_interaction)
                newton_euler_r = sk.cast_NewtonEulerR(inter.relation())
                H = newton_euler_r.jachqT()
                contact_force.append(np.dot(H.T,inter.lambda_(1)))
                u_k=inter.y_k(1)
                u = inter.y(1)
                p = inter.lambda_(1)
                
                work_tot=np.dot(0.5*(u+u_k), p  )
                work_n = 0.5*(u[0]+u_k[0])* p[0]
                work_t = np.dot(0.5*(u+u_k)[1:], p[1:])
                                    
                work.append(work_tot)
              
                if work_n > 2000.:
                    print('u[0]',inter.y(1)[0])
                    print('u_k[0]',inter.y_k(1)[0])
                    print('impact law', u[0]+0.2*u_k[0])
                    print('p[0]',p[0])
                    print('work n ', work_n )
                    #input()

                if work_t > 100. :
                    print('inter ', i_interaction )
                    print('u[0]',inter.y(1)[0])
                    print('u_k[0]',inter.y_k(1)[0])
                    print('impact law', u[0]+0.2*u_k[0])
                    print('p[0]',p[0])
                    idx = np.argwhere(cf_work[:,0] == float(i_interaction))[0][0]
                    print('work n ', work_n, cf_work[idx,1], idx )
                    print('work t ', work_t, cf_work[idx,2])
                    print('cf_work[idx,:]', cf_work[idx,:])
                    idx_cf= np.argwhere(cf[:,22] == float(i_interaction))
                    print('cf', idx_cf, cf[idx_cf,:])
                    #input()

                # if positions is not None:
                #     # print('vold', vold)
                #     # print('v', v)
                #     #print(' u_k = Hvold', np.dot(H,vold))
                #     #print(' u = Hv', np.dot(H,v))
                #     work_new.append(np.dot(0.5*np.dot(H,v+vold), p ))
                #     work_gen= np.dot(0.5*(v+vold),contact_force[-1]) +work_gen
                #idx = idx+1


            #print('work', work)
            # print('work new', work_new)
            # print('work gen', work_gen)
            # print('\n diff cf_work vs. recomputation', np.sum(work)-np.sum(work_new))
            # print('diff cf_work vs. generalized', np.sum(work)- work_gen)

            # p_1 = np.zeros_like(contact_force[0])

            # for f in contact_force:
            #     p_1 = p_1 + f

            # print('p_1', p_1)
            # print('diff P', np.linalg.norm(p_1-neds.p(1)))

            #input()

            normal_work = cf_work[:,1]
            normal_work_negative = np.where( normal_work < 0, normal_work, 0)
            #print('normal_work_negative', normal_work_negative)

            tangent_work = cf_work[:,2]
            tangent_work_negative = np.where( tangent_work < 0, tangent_work, 0)
            #print('tangent_work_negative', tangent_work_negative)


            normal_work_negative_sum = np.sum(normal_work_negative)
            tangent_work_negative_sum = np.sum(tangent_work_negative)
            friction_work_sum = np.sum(cf_work[:,3])

            normal_work_sum = np.sum(normal_work)
            tangent_work_sum = np.sum(tangent_work)
            friction_work_sum = np.sum(cf_work[:,3])



            # print('normal_work_sum', normal_work_sum)
            # print('tangent_work_sum', tangent_work_sum)
            # print('contact_work_sum', normal_work_sum+tangent_work_sum)
            # print('contact_work_sum diff', normal_work_sum+tangent_work_sum-gen_contact_work)
            #print('friction_work_sum', friction_work_sum)

            self._normal_work_sum.append(normal_work_sum)
            self._tangent_work_sum.append(tangent_work_sum)

            self._normal_work_negative_sum.append(normal_work_negative_sum)
            self._tangent_work_negative_sum.append(tangent_work_negative_sum)

            self._friction_work_sum.append(friction_work_sum)
        else:

            self._normal_work_sum.append(0.)
            self._tangent_work_sum.append(0.)

            self._normal_work_negative_sum.append(0.)
            self._tangent_work_negative_sum.append(0.)

            self._friction_work_sum.append(0.)


        normal_work_cum = np.cumsum(self._normal_work_sum)
        tangent_work_cum = np.cumsum(self._tangent_work_sum)
        #print( 'normal_work_cum {0:3e} '.format(normal_work_cum[-1]) )
        #print( 'tangent_work_cum {0:3e}'. format(tangent_work_cum[-1] ))
        #print( 'contact_work_cum {0:3e}'. format(normal_work_cum[-1]+tangent_work_cum[-1] ))


        total_energy = kinetic_sum - external_forces_work_cum[-1] - normal_work_cum[-1] - tangent_work_cum[-1]
        #print('total_energy {0:3e}'.format( total_energy))

        delta_kinetic_energy = kinetic_sum - kinetic_sum_old
        delta_work = external_forces_work_sum + self._normal_work_sum[-1] + self._tangent_work_sum[-1]

        # print('delta_kinetic_energy', delta_kinetic_energy)
        # print('delta_work', delta_work)
        # print('balance', delta_kinetic_energy-delta_work)

        # print('total_energy - old total energy {0:3e}'.format( total_energy-self._total_energy))


        if ((total_energy - self._total_energy)/total_energy > 0.0000001):
            print('old total_energy {0:3e}'.format( self._total_energy))
            print('############################### The total energy has been increased by {0:.3e} ({1:3.2f}%)'.format(total_energy - self._total_energy,(total_energy - self._total_energy)/total_energy*100))

            # input()
        self._total_energy = total_energy

        # if positions is not None:
        #     if np.linalg.norm(neds.p(1)):
        #         input()


        #input()

before_next_step_hook=kinetic_iteration_hook(obj_projectile_id)


#run_options['before_next_step_iteration_hook'] = before_next_step_hook



# Load and run the simulation
with MechanicsHdf5Runner(mode='r+',io_filename=fn) as io:
    io.run(run_options)

# import os
# if Fremond:
#     directory='energy_results_FremondNSL'
# else:
#     directory='energy_results_NewtonNSL'

# if not os.path.exists(directory):
#     os.makedirs(directory)

# import numpy as np
# np.save(os.path.join(directory,'kinetic_energy.npy'), before_next_step_hook._kinetic_sum)
# np.save(os.path.join(directory,'kinetic_energy_wall.npy'), before_next_step_hook._kinetic_sum_wall)
# np.save(os.path.join(directory,'potential_energy.npy'), before_next_step_hook._potential_sum)



# np.save(os.path.join(directory,'normal_work.npy'), before_next_step_hook._normal_work_sum)
# np.save(os.path.join(directory,'tangent_work.npy'), before_next_step_hook._tangent_work_sum)
# np.save(os.path.join(directory,'friction_work.npy'), before_next_step_hook._friction_work_sum)


# np.save(os.path.join(directory,'normal_work_negative.npy'), before_next_step_hook._normal_work_negative_sum)
# np.save(os.path.join(directory,'tangent_work_negative.npy'), before_next_step_hook._tangent_work_negative_sum)
