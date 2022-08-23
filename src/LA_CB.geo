//SetFactory("OpenCASCADE");

//Mesh.CharacteristicLengthFactor = 0.05;  // Coarse
//Mesh.CharacteristicLengthMin = 500;
//Mesh.CharacteristicLengthMax = 1000;
//Mesh.CharacteristicLengthFromCurvature = 1;

Mesh.StlRemoveDuplicateTriangles = 1;
Mesh.Algorithm    = 2; // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) (Default=2)

Mesh.Algorithm3D    = 4; // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree) (Default=1)
Mesh.Recombine3DAll = 1;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.Smoothing = 1;

// Let's merge an STL mesh that we would like to remesh.
Merge "LA_CB.stl";

// We first classify ("color") the surfaces by splitting the original surface
// along sharp geometrical features. This will create new discrete surfaces,
// curves and points
angle = DefineNumber[40, Min 20, Max 120, Step 1, Name "Parameters/Angle for surface detection" ];

// For complex geometries, patches can be too complex, too elongated or too
// large to be parametrized; setting the following option will force the
// creation of patches that are amenable to reparametrization:
forceParametrizablePatches = DefineNumber[0, Choices{0,1}, Name "Parameters/Create surfaces guaranteed to be parametrizable"];

// For open surfaces include the boundary edges in the classification process.
includeBoundary = 1;

// Force curves to be split on given angle:
curveAngle = 180;

ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches, curveAngle * Pi / 180};

// Create a geometry for all the discrete curves and surfaces in the mesh, by
// computing a parametrization for each one
CreateGeometry;
Coherence Mesh;
RefineMesh;

// Create Ellipsoid
cx = 0; cy = 0; cz = 0;

// middle point
cxt = 12700; cyt = 26950; czt = 31620;

//yaw, pitch, roll angles
alpha = Pi/2.5;
beta = 0;
gamma = -1*Pi/12;

// major and minor axis
ma = 6000; mi = 5180;

//Rotation matrix
oneone = Cos(alpha)*Cos(beta);
onetwo = Cos(alpha)*Sin(beta)*Sin(gamma) - Sin(alpha)*Cos(gamma);
onethree = Cos(alpha)*Sin(beta)*Cos(gamma) + Sin(alpha)*Sin(gamma);
twoone = Sin(alpha)*Cos(beta);
twotwo = Sin(alpha)*Sin(beta)*Sin(gamma) + Cos(alpha)*Cos(gamma);
twothree = Sin(alpha)*Sin(beta)*Cos(gamma) - Cos(alpha)*Sin(gamma);
threeone = -1*Sin(beta);
threetwo = Cos(beta)*Sin(gamma);
threethree = Cos(beta)*Cos(gamma);

//Control point
Point(9) = {cxt, cyt, czt};
Point(14) = {(cx*oneone + (cy + mi)*onetwo + cz*onethree) + cxt, (cx*twoone + (cy + mi)*twotwo + cz*twothree) + cyt, (cx*threeone + (cy + mi)*threetwo + cz*threethree) + czt};
Point(10) = {((cx + ma)*oneone + cy*onetwo + cz*onethree) + cxt, ((cx + ma)*twoone + cy*twotwo + cz*twothree) + cyt, ((cx + ma)*threeone + cy*threetwo + cz*threethree) + czt};
Point(11) = {(cx*oneone + cy*onetwo + (cz + ma)*onethree) + cxt, (cx*twoone + cy*twotwo + (cz + ma)*twothree) + cyt, (cx*threeone + cy*threetwo + (cz + ma)*threethree) + czt};
Point(12) = {((cx - ma)*oneone + cy*onetwo + cz*onethree) + cxt, ((cx - ma)*twoone + cy*twotwo + cz*twothree) + cyt, ((cx - ma)*threeone + cy*threetwo + cz*threethree) + czt};
Point(13) = {(cx*oneone + cy*onetwo + (cz - ma)*onethree) + cxt, (cx*twoone + cy*twotwo + (cz - ma)*twothree) + cyt, (cx*threeone + cy*threetwo + (cz - ma)*threethree) + czt};
Point(15) = {(cx*oneone + (cy - mi)*onetwo + cz*onethree) + cxt, (cx*twoone + (cy - mi)*twotwo + cz*twothree) + cyt, (cx*threeone + (cy - mi)*threetwo + cz*threethree) + czt};

// ellipsoid, lines
Ellipse(13) = {10, 9, 11, 11};
Ellipse(14) = {12, 9, 11, 11};
Ellipse(15) = {12, 9, 13, 13};
Ellipse(16) = {13, 9, 10, 10};

Ellipse(17) = {14, 9, 10, 10};
Ellipse(18) = {14, 9, 11, 11};
Ellipse(19) = {14, 9, 12, 12};
Ellipse(20) = {14, 9, 13, 13};

Ellipse(21) = {15, 9, 10, 10};
Ellipse(22) = {15, 9, 11, 11};
Ellipse(23) = {15, 9, 12, 12};
Ellipse(24) = {15, 9, 13, 13};

// line loops for ellipsoid
Line Loop(101) = {16,-17,20};
Surface(101) = {101};
Line Loop(201) = {15, -20, 19};
Surface(201) = {201};
Line Loop(301) = {15, -24, 23};
Surface(301) = {301};
Line Loop(401) = {16, -21, 24};
Surface(401) = {401};
Line Loop(501) = {13, -18, 17};
Surface(501) = {501};
Line Loop(601) = {14, -18, 19};
Surface(601) = {601};
Line Loop(701) = {13, -22, 21};
Surface(701) = {701};
Line Loop(801) = {14, -22, 23};
Surface(801) = {801};

Transfinite Curve {13, 16, 14, 15} = 20 Using Progression 1;
Transfinite Curve {23, 21, 19, 17} = 20 Using Progression 1;
Transfinite Curve {22, 24, 20, 18} = 20 Using Progression 1;

Transfinite Surface {801} = {11, 12, 15};
Transfinite Surface {701} = {11, 15, 10};
Transfinite Surface {301} = {13, 15, 12};
Transfinite Surface {401} = {13, 10, 15};
Transfinite Surface {101} = {13, 14, 10};
Transfinite Surface {501} = {11, 10, 14};
Transfinite Surface {601} = {11, 14, 12};
Transfinite Surface {201} = {13, 12, 14};

// In batch mode the two steps above can be performed with `gmsh t13.stl
// -reparam 40', which will save `t13.msh' containing the parametrizations, and
// which can thus subsequently be remeshed.
 
// Graded mesh option
h_large = 600;
h_small = 400;
radius = 2.5*ma;
Field[1] = Ball;
Field[1].Radius=radius;
Field[1].Thickness = 700;
Field[1].VIn=h_small;
Field[1].VOut=h_large;
Field[1].XCenter = cxt;
Field[1].YCenter = cyt;
Field[1].ZCenter = czt;

Field[2] = Distance;
Field[2].FieldX = 2;
Field[2].FieldY = 2;
Field[2].FieldZ = 2;
Field[2].CurvesList = {7, 8, 9, 10, 11};
Field[2].PointsList = {1, 2, 3, 4, 5};
Field[2].SurfacesList = {2,3,4,5,6,7};
Field[2].NumPointsPerCurve = 200000;

Field[3] = Threshold;
Field[3].InField = 2;
Field[3].DistMax = 2000;
Field[3].DistMin = 1400;
Field[3].SizeMax = 400;
Field[3].SizeMin = 30;

Field[4] = Min;
Field[4].FieldsList = {1,3};
Background Field = 4;

//Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Create a volume as usual

Physical Surface(1) = {2};
Physical Surface(2) = {3};
Physical Surface(3) = {4};
Physical Surface(4) = {5};
Physical Surface(5) = {6};
Physical Surface(6) = {7};
Physical Surface(7) = {101,201,301,401,501,601,701,801}; // Cryoballoon

Surface Loop(11) = {2,3,4,5,6,7, -801,-701,-601,-501,-401,-301,-201,-101};
Volume(1) = {11};

Physical Volume(1) = {1};

//+
Recursive Delete {
  Point{6}; 
}
Recursive Delete {
  Surface{8}; 
}
//+
Coherence;
