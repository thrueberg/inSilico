// sphere radius
R = 0.65;
// sphere centre
cx = 0.5; cy = 0.5; cz = 0.5;

// mesh size
lc = .05;

// centre
Point(1) = {cx,cy,cz,lc};

// axis intersection points -- equator
Point(2) = {cx+R, cy  , cz, lc};
Point(3) = {cx  , cy+R, cz, lc};
Point(4) = {cx-R, cy  , cz, lc};
Point(5) = {cx  , cy-R, cz, lc};
// poles
Point(6) = { cx, cy, cz-R, lc};
Point(7) = { cx, cy, cz+R, lc};

// equator circle
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// meriadians
Circle(5) =  {3,1,6};
Circle(6) =  {6,1,5};
Circle(7) =  {5,1,7};
Circle(8) =  {7,1,3};
Circle(9) =  {2,1,7};
Circle(10) = {7,1,4};
Circle(11) = {4,1,6};
Circle(12) = {6,1,2};

// loops and surfaces
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
// combine all surfaces
Surface Loop(29) = {28,26,16,14,20,24,22,18};
Volume(30) = {29};

