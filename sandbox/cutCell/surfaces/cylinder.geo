// radius
R  = 0.3;
// length
L  = 2.;
// center point
cx = 0.5; cy = 0.5; cz = 0.5;
// direction
dir = 3;

// mesh size
h = 0.05;

// direction vector, bases
dx = 0.; dy = 0.; dz = 0.;
ex = 0.; ey = 0.; ez = 0.;
fx = 0.; fy = 0.; fz = 0.;

If ( dir == 1 )
    dx = 1; ey = 1; fz = 1;
EndIf

If ( dir == 2 )
    dy = 1; ez = 1; fx = 1;
EndIf

If ( dir == 3 )
    dz = 1; ex = 1; fy = 1;
EndIf

// Axis end points
Point(1) = {cx-L/2*dx, cy-L/2*dy, cz-L/2*dz, h};
Point(2) = {cx+L/2*dx, cy+L/2*dy, cz+L/2*dz, h};

// Begin circle points
Point(11) = {cx-L/2*dx+R*ex, cy-L/2*dy+R*ey, cz-L/2*dz+R*ez, h};
Point(12) = {cx-L/2*dx+R*fx, cy-L/2*dy+R*fy, cz-L/2*dz+R*fz, h};
Point(13) = {cx-L/2*dx-R*ex, cy-L/2*dy-R*ey, cz-L/2*dz-R*ez, h};
Point(14) = {cx-L/2*dx-R*fx, cy-L/2*dy-R*fy, cz-L/2*dz-R*fz, h};

// End circle points
Point(21) = {cx+L/2*dx+R*ex, cy+L/2*dy+R*ey, cz+L/2*dz+R*ez, h};
Point(22) = {cx+L/2*dx+R*fx, cy+L/2*dy+R*fy, cz+L/2*dz+R*fz, h};
Point(23) = {cx+L/2*dx-R*ex, cy+L/2*dy-R*ey, cz+L/2*dz-R*ez, h};
Point(24) = {cx+L/2*dx-R*fx, cy+L/2*dy-R*fy, cz+L/2*dz-R*fz, h};

// Circles
Circle(1) = {11,1,12};
Circle(2) = {12,1,13};
Circle(3) = {13,1,14};
Circle(4) = {14,1,11};

Circle(5) = {21,2,22};
Circle(6) = {22,2,23};
Circle(7) = {23,2,24};
Circle(8) = {24,2,21};

// Mantle lines
Line (11) = {11, 21};
Line (12) = {12, 22};
Line (13) = {13, 23};
Line (14) = {14, 24};

// Line Loops
Line Loop (1) = {-1,-2,-3,-4};
Line Loop (2) = {5,6,7,8};

Line Loop (3) = {1,12,-5,-11};
Line Loop (4) = {2,13,-6,-12};
Line Loop (5) = {3,14,-7,-13};
Line Loop (6) = {4,11,-8,-14};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};
Ruled Surface(5) = {5};
Ruled Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
