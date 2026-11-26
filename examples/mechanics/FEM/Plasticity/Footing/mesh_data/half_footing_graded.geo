// Half of a strip footing problem (symmetry axis x = 0): soil [0, L] x [0, L], footing of
// half-width b on the top edge. Graded mesh: element size h_min at the edge of the footing
// (stress singularity) and along it, growing to h_max away from it (gmsh size fields).
// Same physical names as half_footing_2x2_b02_size01.msh (pytripopfem).
// gmsh -2 -format msh22 half_footing_graded.geo -o half_footing_graded.msh
//   (defaults: 441 T3, gmsh 4.14; e.g. -setnumber h_min 0.01 -setnumber L 4 for finer / larger)
If (!Exists(L))
  L = 2;
EndIf
If (!Exists(b))
  b = 0.2;
EndIf
If (!Exists(h_min))
  h_min = 0.02;
EndIf
If (!Exists(h_max))
  h_max = 0.25;
EndIf

Point(1) = {0, 0, 0, h_max};
Point(2) = {L, 0, 0, h_max};
Point(3) = {L, L, 0, h_max};
Point(4) = {b, L, 0, h_max};
Point(5) = {0, L, 0, h_max};

Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right
Line(3) = {3, 4};  // free surface
Line(4) = {4, 5};  // footing
Line(5) = {5, 1};  // symmetry axis

Physical Line("Dirichlet BC1") = {1};
Physical Line("Dirichlet BC2") = {2};
Physical Line("Dirichlet BC3") = {5};
Physical Line("Forces") = {4};

Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = {6};
Physical Surface("Bulk material") = {7};

// size: h_min at the edge of the footing, ~2 h_min along the footing, h_max beyond
Field[1] = Distance;
Field[1].PointsList = {4};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_min;
Field[2].SizeMax = h_max;
Field[2].DistMin = 0.02 * L;
Field[2].DistMax = 0.6 * L;
Field[3] = Distance;
Field[3].CurvesList = {4};
Field[3].Sampling = 50;
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 2 * h_min;
Field[4].SizeMax = h_max;
Field[4].DistMin = 0.05 * L;
Field[4].DistMax = 0.6 * L;
Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
