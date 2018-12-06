//+
Point(1) = {-1, -1, 0, 0.1};
//+
Point(2) = {1, -1, 0, 0.1};
//+
Point(3) = {1, 0, 0, 0.1};
//+
Point(4) = {0, 0, 0, 0.01};
//+
Point(5) = {0, 1, 0, 0.1};
//+
Point(6) = {-1, 1, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
//+
Plane Surface(1) = {1};
//+
Physical Curve("bdry1") = {1, 2, 3, 4, 5};
//+
Physical Curve("bdry2") = {6};
//+
Physical Surface("omega") = {1};