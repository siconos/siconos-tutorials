// Quarter of a strip with a circular hole (x in [0, 0.5], y in [0, 1], hole of radius 0.25
// centred at the origin), graded mesh: element size h_min near the hole and along the
// ligament (bottom edge, where the plastic zone develops), growing to h_max away from them.
// gmsh -2 -format msh22 -setnumber h_min 0.025 strip_with_circular_hole_graded.geo \
//      -o strip_with_circular_hole_graded.msh   (374 triangles, gmsh 4.14)
If (!Exists(h_min))
  h_min = 0.025;
EndIf
If (!Exists(h_max))
  h_max = 0.1;
EndIf

Point(1) = {0.25, 0, 0, h_max};
Point(2) = {0.5, 0, 0, h_max};
Point(3) = {0.5, 1, 0, h_max};
Point(4) = {0, 1, 0, h_max};
Point(5) = {0, 0.25, 0, h_max};
Point(6) = {0, 0, 0, h_max};

Circle(1) = {5, 6, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};

Physical Line("Prescribed displacement") = {4};
Physical Line("Dirichlet BC1") = {2};
Physical Line("Dirichlet BC2") = {5};

Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = {6};
Physical Surface("Bulk material") = {7};

// size field: distance to the hole and to the ligament
Field[1] = Distance;
Field[1].CurvesList = {1, 2};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_min;
Field[2].SizeMax = h_max;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.4;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
