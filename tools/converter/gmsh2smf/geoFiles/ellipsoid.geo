/* -*- mode: c -*- */

// characteristic element length
lc = .05;

// ellipsoid radius
Rx  = 0.3;
Ry  = 0.4;
Rz  = 0.2;

// ellipsoid centre point
Cx = 0.1;
Cy = 0.1;
Cz = -0.1;

// centre of ellipsoid
Point(1) = {Cx, Cy, Cz, lc};

// points in x-y plane on the equator
Point(2) = { Cx+Rx,    Cy,    Cz, lc};
Point(3) = {    Cx, Cy+Ry,    Cz, lc};
Point(4) = { Cx-Rx,    Cy,    Cz, lc};
Point(5) = {    Cx, Cy-Ry,    Cz, lc};

// south and north poles
Point(6) = {    Cx,    Cy, Cz-Rz, lc};
Point(7) = {    Cx,    Cy, Cz+Rz, lc};

// equator ellipse
Ellipse(1) = {2,1,1,3};
Ellipse(2) = {3,1,1,4};
Ellipse(3) = {4,1,1,5};
Ellipse(4) = {5,1,1,2};

// ellipse in y-z plane
Ellipse(5) = {3,1,1,6};
Ellipse(6) = {6,1,1,5};
Ellipse(7) = {5,1,1,7};
Ellipse(8) = {7,1,1,3};

// ellipse in x-z plane
Ellipse(9)  = {2,1,1,7};
Ellipse(10) = {7,1,1,4};
Ellipse(11) = {4,1,1,6};
Ellipse(12) = {6,1,1,2};

// ellipsoid 'triangles' formed by above quarter circles
Line Loop(13) = {2,8,-10};
Ruled Surface(14) = {13};
Line Loop(15) = {10,3,7};
Ruled Surface(16) = {15};
Line Loop(17) = {-8,-9,1};
Ruled Surface(18) = {17};
Line Loop(19) = {-11,-2,5};
Ruled Surface(20) = {19};
Line Loop(21) = {-5,-12,-1};
Ruled Surface(22) = {21};
Line Loop(23) = {-3,11,6};
Ruled Surface(24) = {23};
Line Loop(25) = {-7,4,9};
Ruled Surface(26) = {25};
Line Loop(27) = {-4,12,-6};
Ruled Surface(28) = {27};

// entire ellipsoid surface composed of triangles
Surface Loop(29) = {28,26,16,14,20,24,22,18};

// volume bounded by surface
Volume(30) = {29};

// try also netgen:
// Mesh.Algorithm3D = 4;
