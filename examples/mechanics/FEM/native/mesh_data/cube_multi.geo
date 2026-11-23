

l = 0;
n_cube=50;
For t In {0:n_cube}
    Printf(" t = %g", t);
    Printf(" l = %g", l);
    Point(4*t+1) = {t, 0., 0, 1.0};
    Point(4*t+2) = {t, 1., 0, 1.0};
    Point(4*t+3) = {t, 1., 1., 1.0};
    Point(4*t+4) = {t, 0.0, 1, 1.0};

    Line(l+1) = {4*t+1, 4*t+2};
    Line(l+2) = {4*t+2, 4*t+3};
    Line(l+3) = {4*t+3, 4*t+4};
    Line(l+4) = {4*t+4, 4*t+1};

    Line Loop(l+5) = {l+1, l+2, l+3, l+4};
    Plane Surface(l+6) = {l+5};

    //Line(l+7) = {4*(t-1)+1, 4*t+1};
    //Line(l+8) = {4*(t-1)+2, 4*t+2};
    //Line(l+9) = {4*(t-1)+3, 4*t+3};
    //Line(l+10) = {4*(t-1)+4, 4*t+4};

    // Line Loop(l+11) = {l+1, l+6, l-, l+5};
    // Plane Surface(l+12) = {l+11};


    // Line Loop(l+13) = {l+1, l+6, l-, l+5};
    // Plane Surface(l+14) = {l+13};

    // Line Loop(l+15) = {l+1, l+6, l-, l+5};
    // Plane Surface(l+16) = {l+11};



    l += 6;
EndFor
Printf("end of first loop l = %g", l);

v=0;
Physical Volume("Bulk Material") = {};
Physical Surface("Dirichlet BC") = {6};
Physical Line("Applied Force") = {6*n_cube+3};

For t In {1:n_cube}
    Printf(" t = %g", t);
    Line(l+1) = {4*(t-1)+1, 4*t+1};
    Line(l+2) = {4*(t-1)+2, 4*t+2};
    Line(l+3) = {4*(t-1)+3, 4*t+3};
    Line(l+4) = {4*(t-1)+4, 4*t+4};

    Line Loop(l+5) = {l+1, 6*t+1, -(l+2), -(6*(t-1)+1)};
    Plane Surface(l+6) = {l+5};


    Line Loop(l+7) = {(l+2), (6*t+2), -(l+3), -(6*(t-1)+2)};
    Plane Surface(l+8) = {l+7};

    Line Loop(l+9) = {(l+3), (6*t+3), -(l+4), -(6*(t-1)+3)};
    Plane Surface(l+10) = {l+9};

    Line Loop(l+11) = {(l+4), (6*t+4), -(l+1), -(6*(t-1)+4)};
    Plane Surface(l+12) = {l+11};

    Surface Loop(l+13) = {(l+12), (l+10), l+8,  l+6   ,   6*t+6, 6*(t-1)+6};
    Volume(l+14) = {l+13};

    Physical Volume("Bulk Material") += {l+14};

    l =l + 14;

EndFor
