# This directory contains some simulations of granular materials with contact, impact and friction.



## Simulation of flows of polyhedra into an hopper chute

+ chute_con_rocas.py : Simulation of polyhedrons (convex hulls) into an hopper
  + chute.py : creation of the geometry of the hopper
  + rocas.py : creation of the randomized polyhedrons
+ chute_con_rocas-MoreauJeanGOSI.py : same simulation with the global formulation of the contact problem
+ chute_con_rocas_y_vibradores.py : Simulation with vibrating systems
  + chute_con_vibrador_bottom.py :  vibrating actuators  on the bottom plane of the hopper
  + chute_con_vibrador_middle.py :  vibrating actuators  on the middle plane of the hopper
  + chute_con_vibrador_rear_up.py :  vibrating actuators  on the rear plane of the hopper
  + chute_con_vibrador_top.py : vibrating actuators  on the top plane of the hopper

## Simulation of  a block fragmentation.
+ tesselation.py : simulation of a block fragmentation.
  + n100-id1.tess tesselation file obtained from the NEPER software
	http://neper.sourceforge.net/index.html
  + read_tess.py :  sinple reader of a NEPER .tess file

## Simulation of spheres into a box.
+ spheres_in_a_box.py
+ mkspheres.py : creation of the sample of spheres using LMGC90
http://www.lmgc.univ-montp2.fr/~dubois/LMGC90/Web/Welcome_%21.html


