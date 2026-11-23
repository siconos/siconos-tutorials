#!/bin/bash
#OAR --project siconos
#OAR -l /nodes=3/cpu=1/core=2,walltime=00:01:00
#OAR -n tp1

set -x
# Set path to siconos exe
EXE=/home/$USER/install/siconos/bin/siconos
EX_FILE=/home/$USER/siconos-tutorials/examples/Mechanics/BouncingBall/BouncingBallTS.cpp

# Nodes list
NODES=$(awk -v ORS=, '{print}' "$OAR_FILE_NODES" | sed 's/,$//')
export NODES
set +x

# Load nix env.
# shellcheck source=/dev/null
source /applis/site/nix.sh
set -x
# Number of cores
nbcores=$(wc -l <"$OAR_NODE_FILE")
# Number of nodes
nbnodes=$(sort -u "$OAR_NODE_FILE" | wc -l)
#Name of the first node
firstnode=$(head -1 "$OAR_NODE_FILE")
#Number of cores allocated on the first node (it is the same on all the nodes)
pernode=$(grep -c "$firstnode\$" "$OAR_NODE_FILE")
echo "cores: $nbcores, nodes: $nbnodes, cores on the first node: $pernode"

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
mpirun -np "$(wc -l <"$OAR_FILE_NODES")" --machinefile "$OAR_NODE_FILE" -mca plm_rsh_agent "oarsh" --prefix "$HOME/.nix-profile" "$EXE" "$EX_FILE"
