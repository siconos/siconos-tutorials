cl__1 = 1;
Point(1) = {.25, 0, 0, 1};
Point(2) = {.5, 0, 0, 1};
Point(3) = {.5, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Point(5) = {.0, 0.25, 0, 1};
Point(6) = {.0, 0., 0, 1};


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
