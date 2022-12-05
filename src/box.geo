//Unit 
M = 25;
// Parameters

L1 = 1; E1 = M*2 + 1;
L2 = 1; E2 = M*4 + 1;
L2 = 1; E3 = M*3 + 1;
L2 = 1; E4 = 0.8*M;

x = 1.311;
y = x + 4.0;

SetFactory("OpenCASCADE");

// Extreme corners
Point(1) = {0, 0, 0, 1.0};
Point(2) = {8, 8, 0, 1.0};
Point(3) = {0, 8, 0, 1.0};
Point(4) = {8, 0, 0, 1.0};
Point(5) = {0, 0, 2, 1.0};
Point(6) = {8, 8, 2, 1.0};
Point(7) = {0, 8, 2, 1.0};
Point(8) = {8, 0, 2, 1.0};

// Plane 1
Point(9) = {x, 3, 0, 1.0};
Point(10) = {x, 3, 2, 1.0};
Point(11) = {y, 3, 0, 1.0};
Point(12) = {y, 3, 2, 1.0};

// Plane 2
Point(13) = {x, 5, 0, 1.0};
Point(14) = {x, 5, 2, 1.0};
Point(15) = {y, 5, 0, 1.0};
Point(16) = {y, 5, 2, 1.0};

Line(1) = {16, 15};
Line(2) = {15, 13};
Line(3) = {13, 14};
Line(4) = {14, 16};
Line(5) = {16, 12};
Line(6) = {12, 10};
Line(7) = {10, 9};
Line(8) = {9, 13};
Line(9) = {14, 10};
Line(10) = {11, 15};
Line(11) = {11, 9};
Line(12) = {11, 12};
Line(13) = {2, 3};
Line(14) = {3, 7};
Line(15) = {7, 6};
Line(16) = {6, 2};
Line(17) = {2, 4};
Line(18) = {4, 1};
Line(19) = {1, 5};
Line(20) = {5, 8};
Line(21) = {8, 4};
Line(22) = {8, 6};
Line(23) = {5, 7};
Line(24) = {1, 3};

Curve Loop(1) = {5, 6, -9, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 3, 4, 1};
Plane Surface(2) = {2};
Curve Loop(3) = {8, -2, -10, 11};
Plane Surface(3) = {3};
Curve Loop(4) = {10, -1, 5, -12};
Plane Surface(4) = {4};
Curve Loop(5) = {9, 7, 8, 3};
Plane Surface(5) = {5};
Curve Loop(6) = {6, 7, -11, 12};
Plane Surface(6) = {6};
Curve Loop(7) = {15, -22, -20, 23};
Curve Loop(8) = {4, 5, 6, -9};
Plane Surface(7) = {7, 8};
Curve Loop(9) = {13, -24, -18, -17};
Curve Loop(10) = {2, -8, -11, 10};
Plane Surface(8) = {9, 10};
Curve Loop(11) = {17, -21, 22, 16};
Plane Surface(9) = {11};
Curve Loop(12) = {24, 14, -23, -19};
Plane Surface(10) = {12};
Curve Loop(13) = {21, 18, 19, 20};
Plane Surface(11) = {13};
Curve Loop(14) = {13, 14, 15, 16};
Plane Surface(12) = {14};

Surface Loop(1) = {12, 8, 10, 7, 9, 11, 2, 5, 6, 4};
Volume(1) = {1};
Surface Loop(2) = {2, 5, 6, 4, 3, 1};
Volume(2) = {2};

Transfinite Curve {9, 8, 5, 10, 7, 3, 1, 12} = E1 Using Progression 1;
Transfinite Curve {2, 4, 6, 11} = E2 Using Progression 1;
Transfinite Curve {23, 15, 22, 20, 18, 17, 13, 24} = E3 Using Progression 1;
Transfinite Curve {14, 19, 21, 16} = E4 Using Progression 1;

Transfinite Surface {1} = {16, 14, 10, 12} Right;
Transfinite Surface {3} = {15, 13, 9, 11} Right;
Transfinite Surface {2} = {15, 13, 14, 16};
Transfinite Surface {6} = {11, 9, 10, 12};
Transfinite Surface {4} = {16, 15, 11, 12} Right;
Transfinite Surface {5} = {14, 13, 9, 10} Right;

Transfinite Volume{2} = {11, 12, 10, 9, 15, 16, 14, 13};

Physical Surface("Inlet", 1) = {10};
//+
Physical Surface("Outlet", 2) = {9};
//+
Physical Surface("Bottom", 3) = {11};
//+
Physical Surface("Top", 4) = {12};
//+
Physical Surface("Side1", 5) = {7, 1};
//+
Physical Surface("Side2", 6) = {8, 3};

Physical Volume(1) = {1, 2};

