# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2025 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os
import os.path
from importlib import util
from pathlib import Path

import meshio
import numpy as np

parser = argparse.ArgumentParser(
    description="Postprocess siconos simulation results to generate vtk files"
)
parser.add_argument(
    "--mesh_file",
    help="gmsh file used to build Siconos FE model",
    required=True,
)

parser.add_argument(
    "--simulation_name",
    help="""
    Name of the siconos simulation.
    Usually the name (without ext) of the cpp siconos file
    --> a <name>_displacement.py file must be available in the current dir.
    """,
    required=True,
)

args = parser.parse_args()

mesh_basename = Path(args.mesh_file).stem


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


mesh = meshio.read(
    args.mesh_file,  # string, os.PathLike, or a buffer/open file
    file_format="gmsh",  # optional if filename is a path; inferred from extension
)

# print(mesh.points, mesh.cells, mesh.cells_dict)

# ------------------------------------- vtk output #
# print_meshio_prepost('output mesh in vtk format in ', mesh_basename + ".vtk")
# meshio.write(
#     os.path.splitext(args.mesh_file)[0] + ".vtk",  # str, os.PathLike, or buffer/ open file
#     mesh,
#     # file_format="vtk",  # optional if first argument is a path; inferred from extension
# )

# ------------------------------------- gmsh v2 output #
print_meshio_prepost("output mesh in gmsh v2 format in ", mesh_basename + ".msh2")
meshio.write(
    os.path.splitext(args.mesh_file)[0]
    + ".msh2",  # str, os.PathLike, or buffer/ open file
    mesh,
    file_format="gmsh22",  # optional if first argument is a path; inferred from extension
    binary=False,
    float_fmt=".16e",
)

# ------------------------------------- siconos output #
# print_meshio_prepost('output mesh in cpp format for siconos ', mesh_basename + ".cpp")
# create_Meshcpp_for_siconos(mesh, args.mesh_file)

# ------------------------------------- post processing for paraview"


os.makedirs("vtk", exist_ok=True)


def load_script_as_module(module_name: str, file_path: str):
    """Loads a Python file dynamically as a module."""
    spec = util.spec_from_file_location(module_name, file_path)
    if spec is None:
        raise ImportError(
            f"Could not find spec for module {module_name} at {file_path}"
        )

    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


displacements = load_script_as_module("displacements", "./outputs/T3_displacement.py")
x = displacements.x
y = displacements.y
z = displacements.z

n_samples = len(x)

for i in range(n_samples):
    point_data = {}
    point_data["u_x"] = x[i]
    point_data["u_y"] = y[i]
    point_data["u_z"] = z[i]
    point_data["u"] = np.column_stack((x[i], y[i], z[i]))

    # print('point_data[u]', point_data['u'])

    foutput = "./vtk/" + f"{args.simulation_name}_{i:03d}.vtk"
    print_meshio_prepost("output displacement in vtk format in ", foutput)
    meshio.write_points_cells(
        foutput,
        points=mesh.points,
        cells=mesh.cells,
        # Optionally provide extra data on points, cells, etc.
        point_data=point_data,
        # cell_data=cell_data,
        # field_data=field_data
    )
