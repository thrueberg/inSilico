// Scaffold
 // Cube 1	
	// Points
	Point(1) = {0, 0, 0, 0.01};
	Point(2) = {0, 0, 1, 0.01};
	Point(3) = {0, 1, 0, 0.01};
	Point(4) = {0, 1, 1, 0.01};
	Point(5) = {1, 0, 0, 0.01};
	Point(6) = {1, 0, 1, 0.01};
	Point(7) = {1, 1, 0, 0.01};
	Point(8) = {1, 1, 1, 0.01};

	// Lines
	Line(1) = {1, 3};
	Line(2) = {1, 2};
	Line(3) = {1, 5};
	Line(4) = {2, 4};
	Line(5) = {2, 6};
	Line(6) = {3, 4};
	Line(7) = {3, 7};
	Line(8) = {4, 8};
	Line(9) = {5, 7};
	Line(10) = {5, 6};
	Line(11) = {6, 8};
	Line(12) = {7, 8};

	// Planes
	Line Loop(14) = {6, -4, -2, 1};
	Plane Surface(14) = {14};
	Line Loop(16) = {8, -11, -5, 4};
	Plane Surface(16) = {16};
	Line Loop(18) = {6, 8, -12, -7};
	Plane Surface(18) = {18};
	Line Loop(20) = {12, -11, -10, 9};
	Plane Surface(20) = {20};
	Line Loop(22) = {9, -7, -1, 3};
	Plane Surface(22) = {22};
	Line Loop(24) = {3, 10, -5, -2};
	Plane Surface(24) = {24};

	//Volume
	Surface Loop(26) = {22, 24, 20, 14, 18, 16};
	Volume(26) = {26};

 // Cube 2
	// Points
	Point(21) = {0, 0, 0, 0.01};
	Point(22) = {0, 0, 1, 0.01};
	Point(23) = {0, 1, 0, 0.01};
	Point(24) = {0, 1, 1, 0.01};
	Point(25) = {1, 0, 0, 0.01};
	Point(26) = {1, 0, 1, 0.01};
	Point(27) = {1, 1, 0, 0.01};
	Point(28) = {1, 1, 1, 0.01};

	// Lines
	Line(21) = {21, 23};
	Line(22) = {21, 22};
	Line(23) = {21, 25};
	Line(24) = {22, 24};
	Line(25) = {22, 26};
	Line(26) = {23, 24};
	Line(27) = {23, 27};
	Line(28) = {24, 28};
	Line(29) = {25, 27};
	Line(30) = {25, 26};
	Line(31) = {26, 28};
	Line(32) = {27, 28};

	// Planes
	Line Loop(34) = {26, -24, -22, 21};
	Plane Surface(34) = {14};
	Line Loop(36) = {28, -31, -25, 24};
	Plane Surface(36) = {36};
	Line Loop(38) = {26, 28, -32, -27};
	Plane Surface(38) = {18};
	Line Loop(40) = {32, -31, -30, 29};
	Plane Surface(40) = {40};
	Line Loop(42) = {29, -27, -21, 23};
	Plane Surface(42) = {42};
	Line Loop(44) = {23, 30, -25, -22};
	Plane Surface(44) = {44};

	//Volume
	Surface Loop(46) = {42, 44, 40, 34, 38, 36};
	Volume(46) = {46};


