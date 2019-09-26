#!/bin/bash
#OAR --project siconos
#OAR -l /nodes=3/cpu=1/core=2,walltime=00:01:00
#OAR -n tp1

set -x
# Set path to siconos exe
EXE=/home/$USER/install/siconos/bin/siconos
EX_FILE=/home/$USER/siconos-tutorials/examples/Mechanics/BouncingBall/BouncingBallTS.cpp

# Nodes list
export NODES=`awk -v ORS=, '{print}' $OAR_FILE_NODES|sed 's/,$//'`
set +x

# Load nix env.
source /applis/site/nix.sh
set -x
# Number of cores
nbcores=`cat $OAR_NODE_FILE|wc -l`
# Number of nodes
nbnodes=`cat $OAR_NODE_FILE|sort|uniq|wc -l`
#Name of the first node
firstnode=`head -1 $OAR_NODE_FILE`
#Number of cores allocated on the first node (it is the same on all the nodes)
pernode=`grep "$firstnode\$" $OAR_NODE_FILE|wc -l`
    
####################
# GNU
###################

# -- Nix --
# Nothing to be done. Current profile will be used.

# -- Run mpi --
# options :
#  - nodes list (--machinefile)
#  - send env. var to each node (-x VAR_NAME) 
#  - process number  (-np)
#  - "mca" parameters (see mpi doc), e.g. : -mca <<key>> <<value>> : btl_sm_use_knem ?- 
mpirun -np `cat $OAR_FILE_NODES|wc -l` --machinefile $OAR_NODE_FILE -mca plm_rsh_agent "oarsh" --prefix $HOME/.nix-profile $EXE $EX_FILE


