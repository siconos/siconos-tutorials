"""Helical example,
from B. Caillaud,
http://www.irisa.fr/prive/Benoit.Caillaud/cours-hybride-2014/Cours_modelisation_des_systemes_hybrides/Cours/Cours.html
"""

import numpy as np
import siconos.integrators
import siconos.modeling as sm
import siconos.nonsmooth_formulations
import siconos.numerics as sn
import siconos.plot_config as sicoplot

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)

# == User-defined parameters ==
ndof = 3  # number of degrees of freedom of your system
t0 = 0.0
T = 400  # Total simulation times
h = 0.05  # Time step
rho = 10.0
x00 = rho  # Initial position
x01 = rho
x02 = 0.0
alpha = 0.05  # angle of the square helical
beta = 0.01  # thread of the helical
gamma = 0.0  # thread variation

# -- Dynamical system --
# dx / dt = A.x + b + r
A = np.zeros((ndof, ndof), dtype=np.float64)
x0 = np.zeros(3, dtype=np.float64)
x0.flat[...] = [x00, x01, x02]
b = np.zeros_like(x0)
b[2] = beta
particle = sm.FirstOrderLinearDS(x0, A, b)

# -- Interaction --
# y = C.x + D.lambda
# r = B.lambda
ninter = 2
B = np.zeros((ndof, ninter), dtype=np.float64)
B[0, 0] = -alpha * 0.5
B[0, 1] = 1 + alpha * 0.5
B[1, 0] = -B[0, 1]
B[1, 1] = B[0, 0]
B[2, 0] = gamma * 0.5
B[2, 1] = gamma * 0.5

C = np.zeros((ninter, ndof), dtype=np.float64)
C[0, 0] = C[1, 1] = -1.0

particle_relation = sm.FirstOrderLinearR(C, B)

nslaw = sm.RelayNSL(ninter)

particle_interaction = sm.Interaction(nslaw, particle_relation)

# -- The Model --
filippov = sm.NonSmoothDynamicalSystem(t0, T)
filippov.insertDynamicalSystem(particle)
filippov.link(particle_interaction, particle)

# -- Simulation --
td = siconos.simulation.TimeDiscretisation(t0, h)
simu = siconos.simulation.TimeStepping(filippov, td)
# osi
theta = 0.5
myIntegrator = siconos.integrators.EulerMoreauOSI(theta)
simu.insertIntegrator(myIntegrator)

# osns
osnspb = siconos.nonsmooth_formulations.Relay(sn.solver_ids.SICONOS_RELAY_LEMKE)
simu.insertNonSmoothProblem(osnspb)

# -- Get the values to be plotted --
output_size = 1 + ndof + 2 * ninter
nb_time_steps = int((T - t0) / h) + 1
data_plot = np.empty((nb_time_steps, output_size))
data_plot[0, 0] = filippov.t0()
data_plot[0, 1:4] = particle.x()
data_plot[0, 4:] = 0.0

# time loop
k = 1
while simu.hasNextEvent():
    simu.computeOneStep()
    data_plot[k, 0] = simu.nextTime()
    data_plot[k, 1:4] = particle.x()
    data_plot[k, 4:6] = particle_interaction.lambda_(0)[0:2]
    data_plot[k, 6:8] = particle_interaction.y(0)[0:2]
    k += 1
    simu.nextStep()


# save to disk
np.savetxt("helical.dat", data_plot)


def plot_results():
    """Plot 3d curve z = f(x, y)
    with ds_state = [x, y, z]
    """
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    x = data_plot[:, 1]
    y = data_plot[:, 2]
    z = data_plot[:, 3]
    ax.plot(x, y, z, label="z = f(x, y)")
    ax.legend()
    plt.savefig("helical.png")


# --- Uncomment lines below to plot interesting stuff ---
plot_results()
# plt.subplot(311)
# plt.title('position')
# plt.plot(data_plot[:,0], data_plot[:,1])
# plt.grid()
# plt.subplot(312)
# plt.title('velocity')
# plt.plot(data_plot[:,0], data_plot[:,2])
# plt.grid()
# plt.subplot(313)
# plt.plot(data_plot[:,0], data_plot[:,3])
# plt.title('lambda')
# plt.grid()
# show()

# plt.plot(data_plot[:,1], data_plot[:,2])
# plt.grid()
# show()
