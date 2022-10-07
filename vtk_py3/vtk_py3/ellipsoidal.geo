//Thick-Walled Prolate Spheroid Script

scale = 11;
//lc = 0.05*scale; //this value changes the resolution of the mesh // medium
//lc = 0.025*scale; // Medium 2 mesh
//lc = 0.030*scale; // Medium mesh
lc = <<Meshsize>>*scale; //Coarse
//lz_o = 1*scale; // outer long axis
//lx_o = 0.5*scale; // outer radius
wt=0.15*scale;

lz_o=0.7*0.95*scale;
lx_o=0.35*0.95*scale;
wt=0.112*scale;

lx_i=lx_o-wt; //inner radius
lz_i=lz_o-wt; //inner long axis

Point(1) = {0,0,0,lc};
Point(2) = {lx_o,0,0,lc};
Point(3) = {0,0,-lz_o,lc};
Point(4) = {lx_i,0,0,lc};
Point(5) = {0,0,-lz_i,lc};
Ellipse(1) = {2,1,3,3};
Ellipse(2) = {4,1,5,5};

//OUTER ELLIPSOID

Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{1};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{3};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{6};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{9};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{12};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{15};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{18};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{21};
}

//INNER ELLIPSOID

Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{2};
}

Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{27};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{30};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{33};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{36};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{39};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{42};
}
Extrude {{0,0,1}, {0,0,0}, Pi/4} {
   Line{45};
}

//Upper lines connecting inner and outer half spheroid and the related surfaces


Line(51) = {16, 9};
Line(52) = {15, 8};
Line(53) = {14, 7};
Line(54) = {13, 6};
Line(55) = {4, 2};
Line(56) = {19, 12};
Line(57) = {18, 11};
Line(58) = {17, 10};
Line Loop(59) = {40, 58, -16, -51};
Plane Surface(60) = {59};
Line Loop(61) = {51, -13, -52, 37};
Plane Surface(62) = {61};
Line Loop(63) = {52, -10, -53, 34};
Plane Surface(64) = {63};
Line Loop(65) = {53, -7, -54, 31};
Plane Surface(66) = {65};
Line Loop(67) = {54, -4, -55, 28};
Plane Surface(68) = {67};
Line Loop(69) = {55, -25, -56, 49};
Plane Surface(70) = {69};
Line Loop(71) = {56, -22, -57, 46};
Plane Surface(72) = {71};
Line Loop(73) = {57, -19, -58, 43};
Plane Surface(74) = {73};
Plane Surface(75) = {59};
Surface Loop(76) = {14, 17, 20, 23, 26, 5, 8, 11, 64, 62, 38, 41, 44, 47, 50, 29, 32, 35, 66, 68, 70, 72, 74, 60};
Volume(77) = {76};
Physical Volume(1) = {77};
