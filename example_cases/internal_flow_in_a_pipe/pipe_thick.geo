//+
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0.0, 0.0, 0.0, 3.2, 0, 0, 0.5, 2*Pi};
//+
Cylinder(2) = {0.0, 0.0, 0.0, 3.2, 0, 0, 0.6, 2*Pi};
//+
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }
//+
Transfinite Curve {5, 4} = 60 Using Progression 1;
//+
Transfinite Curve {3, 1} = 80 Using Progression 1;
//+
Transfinite Curve {6, 2} = 80 Using Progression 1;

Physical Surface("Inlet", 1) = {3};
//+
Physical Surface("Outlet", 2) = {2};
//+
Physical Volume("Pipe", 1) = {2};