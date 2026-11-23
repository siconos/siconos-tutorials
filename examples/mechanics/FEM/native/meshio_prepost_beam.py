import os
import os.path
import runpy
from os import path

import meshio

mesh_filename = "./mesh_data/beam.msh"
mesh_filename = "./mesh_data/cube.msh"
mesh_filename = "./mesh_data/cube_multi.msh"
# mesh_filename = './mesh_data/tetra_simple.msh'

mesh_filename = "./mesh_data/square_6.msh"
# mesh_filename = './mesh_data/square_200.msh'
# mesh_filename = './mesh_data/square_2720.msh'
# mesh_filename = './mesh_data/disk.msh'
# mesh_filename = './mesh_data/hertz.msh'

root_mesh_filename = os.path.splitext(mesh_filename)[0]
mesh_basename = os.path.splitext(os.path.split(mesh_filename)[1])[0]


# example_name = 'TH4'
example_name = "beam"
tensor_name = "beamTensor_tensor"

# data written by the simulation (python files defining lists of arrays)
displacement = runpy.run_path(example_name + "_displacement.py")
x, y, z = displacement["x"], displacement["y"], displacement["z"]
tensor = runpy.run_path(tensor_name + ".py")
tensor_xx, tensor_yy, tensor_xy = (
    tensor["tensor_xx"],
    tensor["tensor_yy"],
    tensor["tensor_xy"],
)


def print_meshio_prepost(*args, **kwargs):
    print("[meshio_prepost]", *args, **kwargs)


def create_Meshcpp_for_siconos(mesh, filename):

    with open(os.path.splitext(filename)[0] + ".cpp", "w") as f:
        f.write("#include <Mesh.hpp>\n")
        f.write("static Mesh * createMesh()\n{\n")
        p_cnt = 0
        f.write("std::vector<MeshVertex *> vertices;\n")
        for p in mesh.points:
            f.write(
                f"vertices.push_back(new MeshVertex({p_cnt}, {p[0]}, {p[1]}, {p[2]}));\n"
            )
            p_cnt = p_cnt + 1

        f.write("std::vector<MeshElement *> elements;\n")
        e_cnt = 0
        dim = 0
        for mc in mesh.cells:
            # print("mc.type", mc.type)
            for me in mc[1]:
                # print("me nodes:" ,me)
                str_v = "{"
                for p in me:
                    str_v += f"vertices[{p}],"

                    p_cnt = p_cnt + 1
                str_v += "}"
                # print(str_v)
                f.write(f"std::vector<MeshVertex *> vertices{e_cnt} = {str_v};\n")
                if mc.type == "triangle":
                    type = 2
                    dim = max(dim, 2)
                elif mc.type == "line":
                    type = 1
                    dim = max(dim, 1)
                elif mc.type == "vertex":
                    type = 15
                    dim = max(dim, 0)
                elif mc.type == "tetra":
                    type = 4
                    dim = max(dim, 3)
                f.write(
                    f"elements.push_back(new MeshElement({e_cnt}, {type}, vertices{e_cnt}));\t"
                )

                e_cnt = e_cnt + 1

        f.write(f"Mesh * m =    new Mesh({dim}, vertices, elements);\n")
        f.write("return m;\n }; \n ")


# mesh = meshio.read(
#     mesh_filename,  # string, os.PathLike, or a buffer/open file
#     file_format="gmsh",  # optional if filename is a path; inferred from extension
# )

import numpy as np

# print(mesh.points, mesh.cells, mesh.cells_dict)

# ------------------------------------- vtk output #
# print_meshio_prepost('output mesh in vtk format in ', root_mesh_filename + ".vtk")
# meshio.write(
#     os.path.splitext(mesh_filename)[0] + ".vtk",  # str, os.PathLike, or buffer/ open file
#     mesh,
#     # file_format="vtk",  # optional if first argument is a path; inferred from extension
# )

# ------------------------------------- gmsh v2 output #
# print_meshio_prepost('output mesh in gmsh v2 format in ', root_mesh_filename + ".msh2")
# meshio.write(
#     os.path.splitext(mesh_filename)[0] + ".msh2",  # str, os.PathLike, or buffer/ open file
#     mesh,
#     file_format="gmsh22",  # optional if first argument is a path; inferred from extension
#     binary=False,
#     float_fmt='.16e'
# )

points = [
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 3.0, 0.0],
    [0.0, 4.0, 0.0],
]
cells = [
    ("line", [[0, 1], [1, 2], [2, 3], [3, 4]]),
]

mesh = meshio.Mesh(
    points,
    cells,
)

# ------------------------------------- siconos output #
# print_meshio_prepost('output mesh in cpp format for siconos ', root_mesh_filename + ".cpp")
# create_Meshcpp_for_siconos(mesh, mesh_filename)

# ------------------------------------- post processing for paraview"


if not path.exists("vtk"):
    os.mkdir("vtk")
n_samples = len(x)

for i in range(n_samples):
    point_data = {}
    point_data["u_x"] = x[i]
    point_data["u_y"] = y[i]
    point_data["u_z"] = z[i]
    point_data["u"] = np.column_stack((x[i], y[i], z[i]))
    cell_data = {}
    cell_data["sigmaTension"] = [tensor_xx[i]]
    cell_data["sigmaBendingGP1"] = [tensor_yy[i]]
    cell_data["sigmaBendingGP2"] = [tensor_xy[i]]
    cell_data["sigma"] = [np.column_stack((tensor_xx[i], tensor_yy[i], tensor_xy[i]))]

    print("cell_data[u]", cell_data["sigmaTension"])

    foutput = "./vtk/" + f"{example_name}_{i:03d}.vtk"
    print_meshio_prepost("output displacement in vtk format in ", foutput)
    meshio.write_points_cells(
        foutput,
        points=mesh.points,
        cells=mesh.cells,
        # Optionally provide extra data on points, cells, etc.
        point_data=point_data,
        cell_data=cell_data,
        # field_data=field_data
    )
