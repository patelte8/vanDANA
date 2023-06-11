SetFactory("OpenCASCADE");

// -----------------------------------------------------------

X = 40;

Y1 = 12;
Y2 = 12;
Y3 = 16;

// -----------------------------------------------------------

Cylinder(1) = {0.0, 0.0, 0.0, 1.0, 0, 0, 0.5, 2*Pi};

Cylinder(2) = {0.0, 0.0, 0.0, 1.0, 0, 0, 0.52, 2*Pi};

BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }

// -----------------------------------------------------------

Cylinder(3) = {1, 0, 0, 1, 0, 0, 0.5, 2*Pi};

Cylinder(4) = {1, 0, 0, 1, 0, 0, 0.52, 2*Pi};

BooleanDifference{ Volume{4}; Delete; }{ Volume{3}; Delete; }

// -----------------------------------------------------------

Cylinder(5) = {2.0, 0.0, 0.0, 1.2, 0, 0, 0.5, 2*Pi};

Cylinder(6) = {2.0, 0.0, 0.0, 1.2, 0, 0, 0.52, 2*Pi};

BooleanDifference{ Volume{6}; Delete; }{ Volume{5}; Delete; }

// -----------------------------------------------------------

//+
Coherence;

// -----------------------------------------------------------

Transfinite Curve {5, 4, 3, 1, 9, 7, 13, 11} = X Using Progression 1; 	// circular

Transfinite Curve {6, 2} = Y1 Using Progression 1;

Transfinite Curve {8, 10} = Y2 Using Progression 1;

Transfinite Curve {14, 12} = Y3 Using Progression 1;

// -----------------------------------------------------------

Physical Surface("Inlet", 1) = {3};
Physical Surface("Outlet", 2) = {9};

Physical Volume("Section 1", 1) = {2};
Physical Volume("Section 2", 2) = {4};
Physical Volume("Section 3", 3) = {6};





