//+
Point(1) = {0.0, 0.0, 0, 0.01};
//+
Point(2) = {0.0, 0.5, 0, 1.0};
//+
Point(3) = {0.0, 1., 0, 1.0};
//+
Point(4) = {0.5, 0.5, 0, 0.1};
//+
Point(5) = {-0.5, 0.5, 0, 0.1};
//+
Circle(1) = {5, 2, 3};
//+
Circle(2) = {3, 2, 4};
//+
Circle(3) = {4, 2, 1};
//+
Circle(4) = {1, 2, 5};
//+
Curve Loop(1) = {-1, -2, -3, -4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Applied Force", 5) = {1, 2};
//+
Physical Curve("Contact surface", 4) = {3, 4};
//+
Physical Surface("Bulk Material", 6) = {1};
