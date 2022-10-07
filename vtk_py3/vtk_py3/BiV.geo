Geometry.HideCompounds = 1;

Mesh.CharacteristicLengthFactor = <<Meshsize>> ;// Fine

Mesh.Algorithm    = 6; // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) (Default=2)
Mesh.RecombineAll = 0;

Mesh.RemeshAlgorithm = 1; // (0=no split, 1=automatic, 2=automatic only with metis) (Default=0)

Mesh.RemeshParametrization = 7; // (0=harmonic_circle, 1=conformal_spectral, 2=rbf, 3=harmonic_plane, 4=convex_circle, 5=convex_plane, 6=harmonic square, 7=conformal_fe) (Default=4)

Mesh.Algorithm3D    = 4; // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree) (Default=1)
Mesh.Recombine3DAll = 0;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.Smoothing = 0;

Merge <<LVfilename>>;
Merge <<RVfilename>>;
Merge <<Epifilename>>;

CreateTopology;

ll[] = Line "*";
L_LV_base = newl; Compound Line(L_LV_base) = ll[2];
L_RV_base = newl; Compound Line(L_RV_base) = ll[0];
L_epi_base = newl; Compound Line(L_epi_base) = ll[1];
// Physical Line("EPI_BASE") = {L_epi_base};
//Physical Line(0) = {L_epi_base};

ss[] = Surface "*";
S_LV = news; Compound Surface(S_LV) = ss[0];
S_RV = news; Compound Surface(S_RV) = ss[1];
S_epi = news; Compound Surface(S_epi) = ss[2];
// Physical Surface("LV") = {S_LV};
// Physical Surface("RV") = {S_RV};
// Physical Surface("epi") = {S_epi};
//Physical Surface(3) = {S_epi};
//Physical Surface(1) = {S_RV};
//Physical Surface(2) = {S_LV};


LL_base = newll; Line Loop(LL_base) = {L_LV_base, L_RV_base, L_epi_base};
S_base = news; Plane Surface(S_base) = {LL_base};
// Physical Surface("BASE") = {S_base};
//Physical Surface(4) = {S_base};

SL_wall = newsl; Surface Loop(SL_wall) = {S_LV, S_RV, S_epi, S_base};
V_wall = newv; Volume(V_wall) = {SL_wall};
//Physical Volume("WALL") = {V_wall};
Physical Volume(0) = {V_wall};

