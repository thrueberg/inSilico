//Variables
lc=0.03;
gcl=1;
lcl=0.2;

//Great Cube
	//Points
	Point(101)={0,0,0,lc};
	Point(102)={gcl,0,0,lc};
	Point(103)={0,gcl,0,lc};
	Point(104)={gcl,gcl,0,lc};

	//Lines
	Line(101)={101,102};
	Line(102)={101,103};
	Line(103)={102,104};
	Line(104)={103,104};

	//Surface
	Line Loop(110)={101,103,-104,-102};
	Plane Surface(111)={110};

	//Extrude
	Extrude {0,0,gcl} {Surface{111};}	
	
	Translate {-gcl/2, -gcl/2, -gcl/2} {Volume{1};}
/*
//Lesser Cube
	//Points
	Point(201)={0,0,0,lc};
	Point(202)={lcl,0,0,lc};
	Point(203)={0,lcl,0,lc};
	Point(204)={lcl,lcl,0,lc};

	//Lines
	Line(201)={201,202};
	Line(202)={201,203};
	Line(203)={202,204};
	Line(204)={203,204};

	//Surface
	Line Loop(210)={201,203,-204,-202};
	Plane Surface(211)={210};

	//Extrude
	Extrude {0,0,lcl} {Surface{211};}	
	
	//Translate
	Translate {-lcl/2, -lcl/2, -lcl/2} {Volume{2};}

// Compound volume
	Compound Volume(3) = {1,2};
*/
