

L=12.0;


Point(1) = {L, 0., 0, 1.0};
Point(2) = {L, 1.0, 0, 1.0};
Point(3) = {L, 1.0, 1.0, 1.0};
Point(4) = {L, 0.0, 1.0, 1.0};


Point(5) = {0.0, 0., 0, 1.0};
Point(6) = {0.0, 1.0, 0.0, 1.0};
Point(7) = {0.0, 1.0, 1.0, 1.0};
Point(8) = {0.0, 0.0, 1.0, 1.0};

Point(9) = {L/2., 1.0, 0.0, 1.0};
Point(10) = {L/2., 0.0, 0.0, 1.0};

Point(11) = {L-1., 1.0, 1.0, 1.0};
Point(12) = {L-1., 0.0, 1.0, 1.0};

Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line(5) = {3, 11};
Line(6) = {11, 12};
Line(7) = {12, 4};
Line(8) = {11, 7};
Line(9) = {7, 8};
Line(10) = {8, 12};
Line(11) = {2, 9};
Line(12) = {9, 10};
Line(13) = {10, 1};
Line(14) = {9, 6};
Line(15) = {6, 5};
Line(16) = {5, 10};
Line(17) = {8, 5};
Line(18) = {7, 6};
Line Loop(19) = {13, 3, 11, 12};
Plane Surface(20) = {19};
Line Loop(21) = {3, 4, 1, 2};
Plane Surface(22) = {21};
Line Loop(23) = {7, -1, 5, 6};
Plane Surface(24) = {23};
Line Loop(25) = {10, -6, 8, 9};
Plane Surface(26) = {25};
Line Loop(27) = {18, 15, -17, -9};
Plane Surface(28) = {27};
Line Loop(29) = {14, 15, 16, -12};
Plane Surface(30) = {29};
Line Loop(31) = {8, 18, -14, -11, 4, 5};
Plane Surface(32) = {31};
Line Loop(33) = {10, 7, 2, -13, -16, -17};
Plane Surface(34) = {33};
Surface Loop(35) = {22, 20, 34, 26, 24, 32, 28, 30};
Volume(36) = {35};
Physical Surface("Applied Force") = {24};
Physical Surface("Contact Boundary") = {20};
Physical Surface("Dirichlet BC") = {28};
Physical Volume("Bulk Material") = {36};
