#!/usr/bin/env python

#
# Example of one object under gravity with one contactor and a ground
#

from siconos.mechanics.collision.tools import Contactor
from mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import siconos.mechanics as Mechanics


from numpy import linalg

class Controller(object):

    def __init__(self, io):
        pass

    def initialize(self, io):
        self._io = io
        self._nsds = self._io._nsds
        self._topo = self._nsds.topology()
        self._obj = dict()
        self._nslaw = dict()
        self._interaction = dict()

    def step(self):
        time = self._io._simulation.getTk()
        print('controller - time', time)

        
        nb_ds = self._nsds.getNumberOfDS()
        nb_inter = self._nsds.getNumberOfInteractions()
        
        print('controller - nb_ds', nb_ds)
        print('controller - nb_inter', nb_inter)
        # if (nb_inter > 0):
        #     input()
        list_inter = list()
        for n_inter in range(nb_inter):
            id_inter=self._topo.getInteraction(n_inter).number()
            list_inter.append(id_inter)
            self._interaction[id_inter] = dict()
            self._interaction[id_inter]['interaction object'] = self._topo.getInteraction(n_inter)
            self._interaction[id_inter]['in_contact'] = False

        for k in self._interaction.keys():
            print(k)
            inter= self._interaction[k]['interaction object']
            inter.display()
            _lambda=inter.lambda_(1)
            print(_lambda)
            if (linalg.norm(_lambda) > 1e-15):
                self._interaction[k]['in_contact'] = True


            nsl= Kernel.cast_NewtonImpactFrictionNSL(inter.nonSmoothLaw())
            nsl.display()
            nsl.setMu(0.7)
            
            relation = Kernel.cast_NewtonEulerFrom1DLocalFrameR(inter.relation())
            relation.display()
            print('normal', relation.nc())
            print('pos C1', relation.pc1())
            print('pos C2', relation.pc2())
                #self._obj['DS'] = self._topo.getDynamicalSystem(DS)
      

  





# Creation of the hdf5 file for input/output
with MechanicsHdf5Runner() as io:

    # Definition of a cube as a convex shape
    io.add_convex_shape('Cube', [
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, -1.0, 1.0)])

    # Alternative to the previous convex shape definition.
    # io.add_primitive_shape('Cube1', 'Box', (2, 2, 2))

    # Definition of the ground shape
    io.add_primitive_shape('Ground', 'Box', (100, 100, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.add_Newton_impact_friction_nsl('contact', mu=0.3)

    # The cube object made with an unique Contactor : the cube shape.
    # As a mass is given, it is a dynamic system involved in contact
    # detection and in the simulation.  With no group id specified the
    # Contactor belongs to group 0
    io.add_object('cube', [Contactor('Cube')], translation=[0, 0, 2],
                  velocity=[10, 0, 0, 1, 1, 1],
                  mass=1)

    # the ground object made with the ground shape. As the mass is
    # not given, it is a static object only involved in contact
    # detection.
    io.add_object('ground', [Contactor('Ground')],
                  translation=[0, 0, 0])


# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with MechanicsHdf5Runner(mode='r+') as io:


    C= Controller(io)
 
    
    
    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(with_timer=False,
           gravity_scale=1,
           t0=0,
           T=2.0,
           h=0.0005,
           theta=0.50001,
           Newton_max_iter=20,
           set_external_forces=None,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8,
           numerics_verbose=False,
           output_frequency=None,
           controller =C)
