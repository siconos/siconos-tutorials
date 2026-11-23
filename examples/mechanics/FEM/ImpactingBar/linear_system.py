import time
from math import pi

import numpy as np
import numpy.linalg as LA


class LinearSystem:
    """Description of a mechanical linear sytem with its mass and rigidity matrices and some boundary conditions."""

    def __init__(self, Mass, Stiffness, Bc, fc, fext=None):
        self.Mass = Mass
        self.Stiffness = Stiffness
        ## Boundary conditions : Bc U = fc
        self.Bc = Bc
        self.fc = fc
        self.fext = None
        if fext is not None:
            self.fext = fext

    def computeLinearModes(self):
        """Computation of the eigenfrequencies of a system defined as :
        (K - w^2 M)phi = 0
        """
        # M = L L^t, with L = U sqrt(S) from the singular value decomposition of M
        U, s, _ = LA.svd(self.Mass)
        L = np.dot(U, np.diag(np.sqrt(s)))
        Linv = LA.inv(L)
        # System in its "standard" form, Ax = lambda x, A = L^-1 Stiffness L^-t
        system = np.dot(Linv, np.dot(self.Stiffness, Linv.T))

        t0 = time.perf_counter()
        omeg2, phi = LA.eigh(system)
        print("elapsed with eigh ...", time.perf_counter() - t0)
        freq = np.sqrt(omeg2) / (2.0 * pi)

        # Back to the eigenvectors of the initial system
        vec = np.dot(Linv.T, phi)

        # Check results ...
        # for i in range(vec.shape[1]):
        #    temp = np.dot(Stiffness-omeg2[i]*Mass,vec[:,i])
        #    if(not np.allclose(LA.norm(temp),0.0)):
        #        print "SVD rate ..."

        return freq, vec

    def computeG(self, i, omega):
        """computes G(i,omega) = K - (i*omega)^2 Mass"""

        return self.Stiffness - (i * omega) ** 2 * self.Mass
