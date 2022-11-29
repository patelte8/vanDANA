//Unit 
M = 21;
// Parameters

x = 1.8;
y = x + 1.0;

SetFactory("OpenCASCADE");

// Plane 1
Point(9) = {x, 3.995, 0.5, 1.0};
Point(10) = {x, 3.995, 1.5, 1.0};
Point(11) = {y, 3.995, 0.5, 1.0};
Point(12) = {y, 3.995, 1.5, 1.0};

// Plane 2
Point(13) = {x, 4.005, 0.5, 1.0};
Point(14) = {x, 4.005, 1.5, 1.0};
Point(15) = {y, 4.005, 0.5, 1.0};
Point(16) = {y, 4.005, 1.5, 1.0};

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

Surface Loop(2) = {2, 5, 6, 4, 3, 1};
Volume(1) = {2};

Transfinite Curve {2, 4, 6, 11, 3, 7, 1, 12} = M Using Progression 1;
Transfinite Curve {9, 8, 10, 5} = 7 Using Progression 1;

Transfinite Surface {1} = {16, 14, 10, 12} Right;
Transfinite Surface {3} = {15, 13, 9, 11} Right;
Transfinite Surface {2} = {15, 13, 14, 16};
Transfinite Surface {6} = {11, 9, 10, 12};
Transfinite Surface {4} = {16, 15, 11, 12};
Transfinite Surface {5} = {14, 13, 9, 10};
Transfinite Volume{1} = {11, 12, 10, 9, 15, 16, 14, 13};

Physical Surface("front", 1) = {5};
//+
Physical Volume(1) = {1};
