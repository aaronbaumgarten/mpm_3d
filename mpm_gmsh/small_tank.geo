lc = 0.0025;
Point(1) = {0, 0, 0, lc};
Point(2) = {0.30, 0, 0, lc};
Point(3) = {0.30, 0.1, 0, lc};
Point(4) = {0, 0.1, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {3, 4, 1, 2};
Plane Surface(6) = {6};
Physical Surface(1) = {6};