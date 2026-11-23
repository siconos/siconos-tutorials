import numpy as np
import siconos.plot_config as sicoplot
from mesh import *

# Turn off interactive backend by default
plt, enable_plot = sicoplot.choose_backend(False)
import matplotlib.tri as mtri

# print('coord', coord)
# print('triangle', triangle)

print("number of vertices :", len(coord))
print("number of elements :", len(triangle))

from displacement import *

print("number of displacement samples :", len(x))
print("size of x :", x[0].shape)


coord = np.array(coord)
triangle = np.array(triangle)

plt.figure()
tripost = mtri.Triangulation(coord[:, 0], coord[:, 1], triangles=triangle)
ax = plt.axes()
ax.set_aspect("equal")
ax.triplot(tripost, lw=0.5)
plt.show()

# plt.ylim(-3, 3)
# plt.xlim(0, 3)
# ax= plt.axes()
# ax.set_aspect('equal')
# ax.triplot(tripost,lw = 0.5)

mag = 1e0
k = 0
for x_i in x:
    # print(x[k][:])
    # print(y[k][:])
    tripost = mtri.Triangulation(
        coord[:, 0] + mag * x[k][:], coord[:, 1] + mag * y[k][:], triangles=triangle
    )
    ax = plt.axes()
    # ax.set_aspect('equal')
    ax.triplot(tripost, lw=0.5)
    k = k + 1

plt.figure()
tripost = mtri.Triangulation(coord[:, 0], coord[:, 1], triangles=triangle)
ax = plt.axes()
ax.set_aspect("equal")
ax.triplot(tripost, lw=0.5)
tripost = mtri.Triangulation(
    coord[:, 0] + mag * x[-1][:], coord[:, 1] + mag * y[-1][:], triangles=triangle
)


ax.triplot(tripost, lw=0.5)
echelle_q = ax.tricontourf(tripost, y[-1], cmap=plt.cm.Spectral)
plt.colorbar(echelle_q)
plt.show()

# fig = plt.figure() # initialise la figure
# line, = plt.plot([],[])
# plt.xlim(xmin, xmax)
# plt.ylim(-1,1)


# def init():
#     line.set_data([],[])
#     return line,

# def animate(i):
#     t = i * dt

#     line.set_data(x, y)
#     return line,

# ani = animation.FuncAnimation(fig, animate, init_func=init, frames=100, blit=True, interval=20, repeat=False)
