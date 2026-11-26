Point(1) = {0., 0., 0, 1.0};
Point(2) = {1.0, 0., 0, 1.0};
Point(3) = {0.0, 1.0, 0, 1.0};
Point(4) = {0.0, .0, 1, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};


Line Loop(7) = {5, -4, 1};
Plane Surface(8) = {7};
Line Loop(9) = {1, 2, 3};
Plane Surface(10) = {9};
Line Loop(11) = {6, -4, -3};
Plane Surface(12) = {11};
Line Loop(13) = {6, -5, 2};
Plane Surface(14) = {13};

Surface Loop(15) = {12, 14, 8, 10};
Volume(16) = {15};


Physical Volume("Bulk Material") = {16};
Physical Surface("Dirichlet BC") = {14};
Physical Line("Applied Force") = {1};
