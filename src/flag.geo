// Gmsh project created on Tue Oct 27 14:46:59 2020
SetFactory("OpenCASCADE");

X2 = 10;
X1 = X2*17;
X4 = X2*3;
X5 = X4*2 + X2;

Point(2) = {6.0, 1.9, 0, 1.0};
Point(3) = {6.0, 2.1, 0, 1.0};
Point(4) = {2.489897949, 2.1, 0, 1.0};
Point(5) = {2.489897949, 1.9, 0, 1.0};
Point(6) = {2, 2, 0, 1.0};
Point(7) = {2, 2.5, 0, 1.0};
Point(8) = {2, 1.5, 0, 1.0};

Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {2, 5};

Circle(8) = {4, 6, 7};
Circle(9) = {5, 6, 8};
Circle(10) = {8, 6, 7};
Circle(11) = {5, 6, 4};

Curve Loop(1) = {3, -2, 4, 11};
Plane Surface(1) = {1};
Curve Loop(2) = {8, -10, -9, 11};
Plane Surface(2) = {2};

Transfinite Curve {4} = X1 Using Progression 1;
Transfinite Curve {3} = X1 Using Progression 1;
Transfinite Curve {2} = X2 Using Progression 1;
Transfinite Curve {11} = X2 Using Progression 1;

Transfinite Curve {8, 9} = X4 Using Progression 1;
Transfinite Curve {10} = X5 Using Progression 1;

Transfinite Surface {1};
Coherence;

Physical Curve("cylinder_bd", 2) = {11, 8, 10, 9};
Physical Surface("cylinder", 0) = {2};
Physical Surface("flag", 1) += {1};
