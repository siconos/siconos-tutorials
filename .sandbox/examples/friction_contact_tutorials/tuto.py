# Siconos solvers
import siconos.numerics as sn
# fclib interface
import siconos.fclib as fcl
# h5py
import h5py

import numpy as np
import scipy.linalg as la


# --- Create a friction contact problem ---
# Case 1 : from scratch
# Number of contacts
nc = 3
# W matrix
w_shape = (3 * nc, 3 * nc)
W = np.zeros(w_shape, dtype=np.float64)

W.flat[...] = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]

# q vector
q = np.zeros(3 * nc, dtype=np.float64)
q[...] = [-1, 1, 3, -1, 1, 3, -1, 1, 3]

# Friction coeff
mu = [0.1, 0.1, 0.1]

fc3d = sn.FrictionContactProblem(3, nc, W, q, mu)


# Case 2 : use fclib-library, read hdf5 file




# --- Set solver options ---

# Check Friction_cst.h for a list of solvers ids.
solver_options = sn.SolverOptions(sn.SICONOS_FRICTION_3D_NSGS)#sn.SICONOS_FRICTION_3D_FPP)

eps = np.finfo(np.float64).eps
solver_options.dparam[0] = 100 * eps

# --- Set unknowns ---
velocity = np.zeros(3 * nc, dtype=np.float64)
reaction = np.zeros_like(velocity)

# --- Call solver driver with predefined options ---

sn.fc3d_driver(fc3d, reaction, velocity, solver_options)


