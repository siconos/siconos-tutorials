"""
FEniCS tutorial demo program: Linear elastic problem.

  -div(sigma(u)) = f

The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

from __future__ import print_function
from fenics import *
from ufl import nabla_div
# Scaled variables
L = 1; W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma



print(parameters['linear_algebra_backend'])
print(list_linear_solver_methods())
#parameters['linear_algebra_backend']='cholmod'
input()




# Create mesh and define function space
mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 3, 3, 3)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary condition
tol = 1E-14

def clamped_boundary(x, on_boundary):
    if  (x[0] < tol):
        print('x',x)
        print('on contact boundary')
    return on_boundary and x[0] < tol

bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)

def contact__boundary(x, on_boundary):
    if  (x[2] < tol):
        print('x',x)
        print('on contact boundary')
    return on_boundary and x[2] < tol

contact = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)
print(contact.markers())



# try to play with dofs.

dofmap  = V.dofmap()

print(dofmap.dofs()) # Return list of dof indices on this process that belong to mesh entities of dimension dim
print(dofmap.cell_dofs(0)) # Return list of dof indices on this process that belong to mesh entities of dimension dim
cell = dofmap.cell_dofs(0)

coordinates = V.tabulate_dof_coordinates()
print('coordinates', coordinates)

print('number of nodes : ', coordinates.shape[0])
print('number of dof per  nodes : ', coordinates.shape[1])



input()


# Define strain and stress

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = Function(V)
solve(a == L, u, bc)

#print(u.array())
print(u.compute_vertex_values()) # values order in the numbering of vertex.


vertex_values = u.compute_vertex_values()
for i, x in enumerate(coordinates):
    print('vertex %d: vertex_values[%d]'%(i,i) )
    print('x=',x)
    print('u(x)', u(x))
    #print('vertex %d: vertex_values[%d] = %g\tu(%s) = %g' %
      #    (i, i, vertex_values[i], x, u(x)))

element = V.element()
dofmap = V.dofmap()
      
for cell in cells(mesh):
    print('cell number=',  cell.index())
    
    print(element.tabulate_dof_coordinates(cell))
    print(dofmap.cell_dofs(cell.index()))

      
input()


# Plot solution
plot(u, title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plot(von_Mises, title='Stress intensity')

print(u.vector().array_view())

F = assemble(L)
print(F.get_local())

K_mat= assemble(a)
print(K_mat.array())

import numpy as np

K_mat_np = np.array(K_mat.array())
eig, eig_v = np.linalg.eig(K_mat_np)
print('eigenvalues before bc', eig.min(), eig.max())
input()

bc.apply(K_mat)
K_mat_np = np.array(K_mat.array())
eig, eig_v = np.linalg.eig(K_mat_np)
print('eigenvalues after bc', eig.min(), eig.max())
input()

bc.apply(F)

input()


print(as_backend_type(K_mat))


#print(K_mat.getValuesCSR()[::-1])

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
#plot(u_magnitude, 'Displacement magnitude')
print('min/max u:',
       u_magnitude.vector().array_view().min(),
       u_magnitude.vector().array_view().max())



# Save solution to file in VTK format
File('elasticity/displacement.pvd') << u
File('elasticity/von_mises.pvd') << von_Mises
File('elasticity/magnitude.pvd') << u_magnitude

# Hold plot
#interactive()
