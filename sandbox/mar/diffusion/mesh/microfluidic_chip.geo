/*********************************************************************
 *
 *  Gmsh Microfluidic device
 *
 *  Chip geometry parameterization with characteristic lengths variable 
 *
 *********************************************************************/
/////////////////////////////////////////
//             INPUT VAR              //
////////////////////////////////////////
//
// See attached figure on    .jpg
//
// ncol1=number of posts(columns) on each side of the gel region
// hcol1=height of the posts
// hgel1=height of the gel region
// bgel1=width of the gel region
// bcan1=width of the control and condition channel (assumption of equal width)
// rdep1=input reservoir radius 
// rdep2=output reservoir radius 
//
// The variables values are on mm
// But you can also introduce on microns, but all on the same units
//
ncol1=6;
hcol1=0.29;
hgel1=3.22;

bgel1=1.3;

bcan1=0.92;
rdep1=2.25;
rdep2=2.25;
//
// Characteristic element length 
//
cl1 = 0.05;

////////////////////////////////
//          3D MODEL          //
////////////////////////////////
//
// chip depth (if required)
//
pchip1=0.146;


//////////////////////////////////////
//        POINTS DEFINITION        //
/////////////////////////////////////
//
// 1.- In this section are defined the points that generate the outline of the gel
//
// p1 represents the origin of the coordinate system 
// p1-p4 are the points corresponding with the upper left part of the gel contour
p1 = newp; Point(p1) = {0, 0, 0, cl1};
p2 = newp; Point(p2) = {-bgel1*3/13, -hgel1*9/161, 0, cl1};
p3 = newp; Point(p3) = {-bgel1*3/13, -hgel1*12/161, 0, cl1};
p4x = -bgel1*5/13;
p4y = -hgel1*18/161;
p4 = newp; Point(p4) = {p4x, p4y , 0, cl1};

// p25-p28 are the points corresponding with the bottom left part of the gel contour
p25x = -bgel1*5/13;
p25y = -hgel1*143/161;
p25 = newp; Point(p25) = {p25x, p25y, 0, cl1};
p26 = newp; Point(p26) = {-bgel1*3/13, -hgel1*149/161, 0, cl1};
p27 = newp; Point(p27) = {-bgel1*3/13, -hgel1*152/161, 0, cl1};
p28 = newp; Point(p28) = {0, -hgel1, 0, cl1};

// p29-p32 are the points corresponding with the upper right part of the gel contour
p29 = newp; Point(p29) = {bgel1*3/13, 0, 0, cl1};
p30 = newp; Point(p30) = {bgel1*6/13, -hgel1*9/161, 0, cl1};
p31 = newp; Point(p31) = {bgel1*6/13, -hgel1*12/161, 0, cl1};
p32x = bgel1*8/13;
p32y = -hgel1*18/161;
p32 = newp; Point(p32) = {p32x, p32y, 0, cl1};

// p53-p56 are the points corresponding with the bottom right part of the gel contour
p53x = bgel1*8/13;
p53y = -hgel1*143/161;
p53 = newp; Point(p53) = {p53x, p53y, 0, cl1};
p54 = newp; Point(p54) = {bgel1*6/13, -hgel1*149/161, 0, cl1};
p55 = newp; Point(p55) = {bgel1*6/13, -hgel1*152/161, 0, cl1};
p56 = newp; Point(p56) = {bgel1*3/13, -hgel1, 0, cl1};




//////////////////////////////////////////////////////////////////////////////////////////
// In the following sections 2, 3, 4 and 5 have been defined the intersection of the 
// channels with the reservoirs by the intersection equation between a line and a circle. 
//
// Obtaining for each case six parameters a, b​​, c, d, e and f in terms of the variables 
// of the problem, where m is the slope at the end of the channel.
/////////////////////////////////////////////////////////////////////////////////////////
//
// 2.- In this section are defined the points that generate the outline of the upper part  
//    of the control channel
//
p100 = newp; Point(p100) = {p4x-bcan1, p4y, 0, cl1};
p101 = newp; Point(p101) = {p4x-227/200*bcan1, p4y, 0, cl1};
p102 = newp; Point(p102) = {p4x-2399/2300*bcan1, p4y+91/920*bcan1, 0, cl1};
p103 = newp; Point(p103) = {p4x-2219/920*bcan1, p4y+6327/4600*bcan1, 0, cl1};
p104 = newp; Point(p104) = {p4x-22909/9200*bcan1, p4y+12887/9200*bcan1, 0, cl1};
p105x = p4x-70/23*bcan1;
p105y = p4y+489/368*bcan1;
p105 = newp; Point(p105) = {p105x, p105y, 0, cl1};
p106 = newp; Point(p106) = {p4x-3333/9200*bcan1, p4y+478/575*bcan1, 0, cl1};
p107 = newp; Point(p107) = {p4x-15927/9200*bcan1, p4y+1212/575*bcan1, 0, cl1};
p108 = newp; Point(p108) = {p4x-2403/920*bcan1, p4y+11009/4600*bcan1, 0, cl1};
p109x = p4x-29421/9200*bcan1;
p109y = p4y+21357/9200*bcan1;
p109 = newp; Point(p109) = {p109x, p109y, 0, cl1};
p112 = newp; Point(p112) = {p4x-1443/460*bcan1, p4y+16791/9200*bcan1, 0, cl1};
p113 = newp; Point(p113) = {p4x-57/23*bcan1, p4y+30/23*bcan1, 0, cl1};
p114x = p4x-1443/460*bcan1-rdep1*Cos( 7*Pi/180 );
p114y = p4y+16791/9200*bcan1-rdep1*Sin( 7*Pi/180 );
p114 = newp; Point(p114) = {p114x, p114y, 0, cl1};
p115x = p114x-rdep1*Cos( 7*Pi/180 );
p115y = p114y-rdep1*Sin( 7*Pi/180 );
p115 = newp; Point(p115) = {p115x, p115y, 0, cl1};
//
// Calculation of the control channel cutoffs with the upper input reservoir: 
// p110 and p111
//
	m1 = 661/5391;
	e1 = p109y-m1*p109x;
	f1 = e1-p114y;
//
	a1 = 1+(m1^2);
	b1 = (2*m1*f1)-(2*p114x);
	c1 = (p114x^2)+(f1^2)-(rdep1^2);
	p110x = (-b1+Sqrt(b1^2-(4*a1*c1)))/(2*a1);
	p110y = (m1*p110x)+e1;
	p110 = newp; Point(p110) = {p110x, p110y, 0, cl1};
//
	e2 = p105y-m1*p105x;
	f2 = e2-p114y;
//
	a2 = 1+(m1^2);
	b2 = (2*m1*f2)-(2*p114x);
	c2 = (p114x^2)+(f2^2)-(rdep1^2);
	p111x = (-b2+Sqrt((b2^2)-(4*a2*c2)))/(2*a2);
	p111y = (m1*p111x)+e2;
	p111 = newp; Point(p111) = {p111x, p111y, 0, cl1};
//
// 3.- In this section are defined the points that generate the outline of the bottom part  
//    of the control channel
//
	p200 = newp; Point(p200) = {p25x-bcan1, p25y, 0, cl1};
	p201 = newp; Point(p201) = {p25x-227/200*bcan1, p25y, 0, cl1};
	p202 = newp; Point(p202) = {p25x-2399/2300*bcan1, p25y-91/920*bcan1, 0, cl1};
	p203 = newp; Point(p203) = {p25x-2219/920*bcan1, p25y-6327/4600*bcan1, 0, cl1};
	p204 = newp; Point(p204) = {p25x-22909/9200*bcan1, p25y-12887/9200*bcan1, 0, cl1};
	p205x = p25x-70/23*bcan1;
	p205y = p25y-489/368*bcan1;
	p205 = newp; Point(p205) = {p205x, p205y, 0, cl1};
	p206 = newp; Point(p206) = {p25x-3333/9200*bcan1, p25y-478/575*bcan1, 0, cl1};
	p207 = newp; Point(p207) = {p25x-15927/9200*bcan1, p25y-1212/575*bcan1, 0, cl1};
	p208 = newp; Point(p208) = {p25x-2403/920*bcan1, p25y-11009/4600*bcan1, 0, cl1};
	p209x = p25x-29421/9200*bcan1;
	p209y = p25y-21357/9200*bcan1;
	p209 = newp; Point(p209) = {p209x, p209y, 0, cl1};
	p212 = newp; Point(p212) = {p25x-1443/460*bcan1, p25y-16791/9200*bcan1, 0, cl1};
	p213 = newp; Point(p213) = {p25x-57/23*bcan1, p25y-30/23*bcan1, 0, cl1};
	p214x = p25x-1443/460*bcan1-rdep1*Cos( 7*Pi/180 );
	p214y = p25y-16791/9200*bcan1+rdep1*Sin( 7*Pi/180 );
	p214 = newp; Point(p214) = {p214x, p214y, 0, cl1};
	p215x = p214x-rdep1*Cos( 7*Pi/180 );
	p215y = p214y+rdep1*Sin( 7*Pi/180 );
	p215 = newp; Point(p215) = {p215x, p215y, 0, cl1};
//
// Calculation of the control channel cutoffs with the lower input reservoir: 
// p210 and p211
//
	m3 = -661/5391;
	e3 = p209y-m3*p209x;
	f3 = e3-p214y;
//
	a3 = 1+(m3^2);
	b3 = (2*m3*f3)-(2*p214x);
	c3 = (p214x^2)+(f3^2)-(rdep1^2);
	p210x = (-b3+Sqrt(b3^2-(4*a3*c3)))/(2*a3);
	p210y = (m3*p210x)+e3;
	p210 = newp; Point(p210) = {p210x, p210y, 0, cl1};
//
	e4 = p205y-m3*p205x;
	f4 = e4-p214y;
//
	a4 = 1+(m3^2);
	b4 = (2*m3*f4)-(2*p214x);
	c4 = (p214x^2)+(f4^2)-(rdep1^2);
	p211x = (-b4+Sqrt((b4^2)-(4*a4*c4)))/(2*a4);
	p211y = (m3*p211x)+e4;
	p211 = newp; Point(p211) = {p211x, p211y, 0, cl1};
//
//
// 4.- In this section are defined the points that generate the outline of the upper part  
//    of the condition channel
//
	p300 = newp; Point(p300) = {p32x+bcan1, p32y, 0, cl1};
	p301 = newp; Point(p301) = {p32x+227/200*bcan1, p32y, 0, cl1};
	p302 = newp; Point(p302) = {p32x+2399/2300*bcan1, p32y+91/920*bcan1, 0, cl1};
	p303 = newp; Point(p303) = {p32x+2219/920*bcan1, p32y+6327/4600*bcan1, 0, cl1};
	p304 = newp; Point(p304) = {p32x+22909/9200*bcan1, p32y+12887/9200*bcan1, 0, cl1};
	p305x = p32x+70/23*bcan1;
	p305y = p32y+489/368*bcan1;
	p305 = newp; Point(p305) = {p305x, p305y, 0, cl1};
	p306 = newp; Point(p306) = {p32x+3333/9200*bcan1, p32y+478/575*bcan1, 0, cl1};
	p307 = newp; Point(p307) = {p32x+15927/9200*bcan1, p32y+1212/575*bcan1, 0, cl1};
	p308 = newp; Point(p308) = {p32x+2403/920*bcan1, p32y+11009/4600*bcan1, 0, cl1};
	p309x = p32x+29421/9200*bcan1;
	p309y = p32y+21357/9200*bcan1;
	p309 = newp; Point(p309) = {p309x, p309y, 0, cl1};
	p312 = newp; Point(p312) = {p32x+1443/460*bcan1, p32y+16791/9200*bcan1, 0, cl1};
	p313 = newp; Point(p313) = {p32x+57/23*bcan1, p32y+30/23*bcan1, 0, cl1};
	p314x = p32x+1443/460*bcan1+rdep2*Cos( 7*Pi/180 );
	p314y = p32y+16791/9200*bcan1-rdep2*Sin( 7*Pi/180 );
	p314 = newp; Point(p314) = {p314x, p314y, 0, cl1};
	p315x = p314x+rdep2*Cos( 7*Pi/180 );
	p315y = p314y-rdep2*Sin( 7*Pi/180 );
	p315 = newp; Point(p315) = {p315x, p315y, 0, cl1};
//
// Calculation of the control channel cutoffs with the upper output reservoir: 
// p310 and p311
//
	e5 = p309y-m3*p309x;
	f5 = e5-p314y;
//
	a5 = 1+(m3^2);
	b5 = (2*m3*f5)-(2*p314x);
	c5 = (p314x^2)+(f5^2)-(rdep2^2);
	p310x = (-b5-Sqrt(b5^2-(4*a5*c5)))/(2*a5);
	p310y = (m3*p310x)+e5;
	p310 = newp; Point(p310) = {p310x, p310y, 0, cl1};
//
	e6 = p305y-m3*p305x;
	f6 = e6-p314y;
//
	a6 = 1+(m3^2);
	b6 = (2*m3*f6)-(2*p314x);
	c6 = (p314x^2)+(f6^2)-(rdep2^2);
	p311x = (-b6-Sqrt((b6^2)-(4*a6*c6)))/(2*a6);
	p311y = (m3*p311x)+e6;
	p311 = newp; Point(p311) = {p311x, p311y, 0, cl1};
//
// 5.- In this section are defined the points that generate the outline of the bottom part  
//    of the condition channel
//
	p400 = newp; Point(p400) = {p53x+bcan1, p53y, 0, cl1};
	p401 = newp; Point(p401) = {p53x+227/200*bcan1, p53y, 0, cl1};
	p402 = newp; Point(p402) = {p53x+2399/2300*bcan1, p53y-91/920*bcan1, 0, cl1};
	p403 = newp; Point(p403) = {p53x+2219/920*bcan1, p53y-6327/4600*bcan1, 0, cl1};
	p404 = newp; Point(p404) = {p53x+22909/9200*bcan1, p53y-12887/9200*bcan1, 0, cl1};
	p405x = p53x+70/23*bcan1;
	p405y = p53y-489/368*bcan1;
	p405 = newp; Point(p405) = {p405x, p405y, 0, cl1};
	p406 = newp; Point(p406) = {p53x+3333/9200*bcan1, p53y-478/575*bcan1, 0, cl1};
	p407 = newp; Point(p407) = {p53x+15927/9200*bcan1, p53y-1212/575*bcan1, 0, cl1};
	p408 = newp; Point(p408) = {p53x+2403/920*bcan1, p53y-11009/4600*bcan1, 0, cl1};
	p409x = p53x+29421/9200*bcan1;
	p409y = p53y-21357/9200*bcan1;
	p409 = newp; Point(p409) = {p409x, p409y, 0, cl1};
	p412 = newp; Point(p412) = {p53x+1443/460*bcan1, p53y-16791/9200*bcan1, 0, cl1};
	p413 = newp; Point(p413) = {p53x+57/23*bcan1, p53y-30/23*bcan1, 0, cl1};
	p414x =p53x+1443/460*bcan1+rdep2*Cos( 7*Pi/180 );
	p414y =-hgel1*143/161-16791/9200*bcan1+rdep2*Sin( 7*Pi/180 );
	p414 = newp; Point(p414) = {p414x, p414y, 0, cl1};
	p415x = p414x+rdep2*Cos( 7*Pi/180 );
	p415y = p414y+rdep2*Sin( 7*Pi/180 );
	p415 = newp; Point(p415) = {p415x, p415y, 0, cl1};
//
// Calculation of the control channel cutoffs with the lower output reservoir: 
// p410 and p411
//
	e7 = p409y-m1*p409x;
	f7 = e7-p414y;
//
	a7 = 1+(m1^2);
	b7 = (2*m1*f7)-(2*p414x);
	c7 = (p414x^2)+(f7^2)-(rdep2^2);
	p410x = (-b7-Sqrt(b7^2-(4*a7*c7)))/(2*a7);
	p410y = (m1*p410x)+e7;
	p410 = newp; Point(p410) = {p410x, p410y, 0, cl1};
//
	e8 = p405y-m1*p405x;
	f8 = e8-p414y;
//
	a8 = 1+(m1^2);
	b8 = (2*m1*f8)-(2*p414x);
	c8 = (p414x^2)+(f8^2)-(rdep2^2);
	p411x = (-b8-Sqrt((b8^2)-(4*a8*c8)))/(2*a8);
	p411y = (m1*p411x)+e8;
	p411 = newp; Point(p411) = {p411x, p411y, 0, cl1};
//
//
//////////////////////////////////////
//         LINES DEFINITION        //
/////////////////////////////////////
//
//
// In this section are defined the lines that generate the outline of the Microfluidic device.
//
//
//////////////////////////////////////////////////////////////////////////////
// Condition that controls that input and output reservoirs do not overlap
// In case that the overlap exists does not draw the microfluidic device
// and exits.
//////////////////////////////////////////////////////////////////////////////
//
    If (Fabs(p114y)+Fabs(p214y) >= (rdep1+rdep2+0.05))
//
// The contour lines of the gel are generated from the points defined above
//
		Line(1) = {p1, p2};
		Line(2) = {p2, p3};
		Line(3) = {p3, p4};
		Line(19) = {p25, p26};
		Line(20) = {p26, p27};
		Line(21) = {p27, p28};
		Line(22) = {p1, p29};
		Line(23) = {p29, p30};
		Line(24) = {p30, p31};
		Line(25) = {p31, p32};
		Line(41) = {p53, p54};
		Line(42) = {p54, p55};
		Line(43) = {p55, p56};
		Line(44) = {p56, p28};
		Line(55) = {p100, p200};
		Line(56) = {p300, p400};
//
// Now the contour lines of the top left side of the microfluidic device are generated
// from the points defined above
//
		Line(100) = {p102, p103};
		Line(101) = {p104, p105};
		Line(102) = {p105, p111};
		Line(103) = {p106, p107};
		Line(104) = {p108, p109};
		Line(105) = {p109, p110};
		Circle(1000) = {p100, p101, p102};
		Circle(1001) = {p4, p101, p106};
		Circle(1002) = {p103, p113, p104};
		Circle(1003) = {p107, p113, p108};
		Circle(1004) = {p115, p114, p110};
		Circle(1005) = {p115, p114, p111};
//
// Now the contour lines of the bottom left side of the microfluidic device are generated
// from the points defined above
//
		Line(200) = {p202, p203};
		Line(201) = {p204, p205};
		Line(202) = {p205, p211};
		Line(203) = {p206, p207};
		Line(204) = {p208, p209};
		Line(205) = {p209, p210};
		Circle(2000) = {p200, p201, p202};
		Circle(2001) = {p25, p201, p206};
		Circle(2002) = {p203, p213, p204};
		Circle(2003) = {p207, p213, p208};
		Circle(2004) = {p215, p214, p210}; 
		Circle(2005) = {p215, p214, p211};
//
// Now the contour lines of the top right side of the microfluidic device are generated
// from the points defined above
//
		Line(300) = {p302, p303};
		Line(301) = {p304, p305};
		Line(302) = {p305, p311};
		Line(303) = {p306, p307};
		Line(304) = {p308, p309};
		Line(305) = {p309, p310};
		Circle(3000) = {p300, p301, p302};
		Circle(3001) = {p32, p301, p306};
		Circle(3002) = {p303, p313, p304};
		Circle(3003) = {p307, p313, p308};
		Circle(3004) = {p315, p314, p310};
		Circle(3005) = {p315, p314, p311};
//
// Now the contour lines of the bottom right side of the microfluidic device are generated
// from the points defined above
//
		Line(400) = {p402, p403};
		Line(401) = {p404, p405};
		Line(402) = {p405, p411};
		Line(403) = {p406, p407};
		Line(404) = {p408, p409};
		Line(405) = {p409, p410};
		Circle(4000) = {p400, p401, p402};
		Circle(4001) = {p53, p401, p406};
		Circle(4002) = {p403, p413, p404};
		Circle(4003) = {p407, p413, p408};
		Circle(4004) = {p415, p414, p410};
		Circle(4005) = {p415, p414, p411};
//
//
//
//////////////////////////////////////
//         POSTS DEFINITION        //
/////////////////////////////////////
//
//
// First of all, the calculation of these three parameters, that will define the 
// position of the posts on the gel region as well as their separation.
//
// hdif1=the total height corresponding to the contact area between the channel and the gel,
//       this zone is the total height where the substance may diffuse.
// hcoltot1=the height corresponding to one post
// hnocol1=the height corresponding to the separation between two posts
//
//
		hdif1 = 125/161*hgel1;
		hcoltot1 = hcol1*ncol1;
		hnocol1 = (hdif1-hcoltot1)/(ncol1+1);
//
// FIRST CONDITION: (not included)
// 		altura mínima que debe haber entre las columnas para que las células
// 		sean capaces de pasar al gel ha de ser hnocol=0.05 mm=50 micras
//
//
// SECOND CONDITION: (included)
//  	the total height contribution by the posts, user-defined, is less than the section defined for it.
//  	We defined the minimum separation between each post equal to 1 micron.
//		minhnocoltot1=minimum total height between posts 
//
		minhnocoltot1=(ncol1+1)*0.01;

		If ((hcoltot1+minhnocoltot1) <= hdif1)
//
// Loop for that creates the two columns of posts and the gaps between them depending on the number of 
// posts and the height selected. The width and the output slope of the posts is fixed and it is assumed 
// equal for all of them.
//
			contador=0;
			t = ncol1;

			For t In {1:ncol1}

// 1.- In this section are defined the points that generate the posts
			
				If (t==1)
					p1000 = newp; Point(p1000) = {p4x, p4y-(hnocol1*t)-(hcol1*(t-1)), 0, cl1};
					p1001 = newp; Point(p1001) = {p4x, p4y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p1002 = newp; Point(p1002) = {p4x, p4y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p1003 = newp; Point(p1003) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1};
					p1004 = newp; Point(p1004) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*t)+(16/644*hgel1), 0, cl1};
		
					p2000 = newp; Point(p2000) = {p32x, p32y-(hnocol1*t)-(hcol1*(t-1)), 0, cl1};
					p2001 = newp; Point(p2001) = {p32x, p32y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p2002 = newp; Point(p2002) = {p32x, p32y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p2003 = newp; Point(p2003) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1} ;
					p2004 = newp; Point(p2004) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*t)+16/644*hgel1, 0, cl1};
				EndIf
	
				If (t>1 && t<ncol1)
					p1000 = p1002;
					p1001 = newp; Point(p1001) = {p4x, p4y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p1002 = newp; Point(p1002) = {p4x, p4y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p1003 = newp; Point(p1003) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1};
					p1004 = newp; Point(p1004) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*t)+(16/644*hgel1), 0, cl1};
		
					p2000 = p2002;
					p2001 = newp; Point(p2001) = {p32x, p32y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p2002 = newp; Point(p2002) = {p32x, p32y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p2003 = newp; Point(p2003) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1} ;
					p2004 = newp; Point(p2004) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*t)+16/644*hgel1, 0, cl1};
				EndIf
	
				If (t==ncol1)
					p1000 = p1002;
					p1001 = newp; Point(p1001) = {p4x, p4y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p1002 = newp; Point(p1002) = {p4x, p4y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p1003 = newp; Point(p1003) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1};
					p1004 = newp; Point(p1004) = {-bgel1*18/65, p4y-(hnocol1*t)-(hcol1*t)+(16/644*hgel1), 0, cl1};
		
					p2000 = p2002;
					p2001 = newp; Point(p2001) = {p32x, p32y-(hnocol1*t)-(hcol1*t), 0, cl1};
					p2002 = newp; Point(p2002) = {p32x, p32y-(hnocol1*(t+1))-(hcol1*t), 0, cl1};
					p2003 = newp; Point(p2003) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*(t-1))-(16/644*hgel1), 0, cl1} ;
					p2004 = newp; Point(p2004) = {bgel1*33/65, p32y-(hnocol1*t)-(hcol1*t)+16/644*hgel1, 0, cl1};
				EndIf
				
// 2.- In this section are defined the lines that generate the posts
    
				If (t==1)
					Line(57) = {p4, p1000};
					Line(58) = {p32, p2000};
					l1 = newl; Line(l1) = {p1001, p1000};
					l2 = newl; Line(l2) = {p1001, p1002};
					l3 = newl; Line(l3) = {p2001, p2000}; 
					l4 = newl; Line(l4) = {p2001, p2002}; 
			
					l5 = newl; Line(l5) = {p1000, p1003};
					l6 = newl; Line(l6) = {p1003, p1004};
					l7 = newl; Line(l7) = {p1004, p1001}; 
					l8 = newl; Line(l8) = {p2000, p2003}; 
					l9 = newl; Line(l9) = {p2003, p2004}; 
					l10 = newl; Line(l10) = {p2004, p2001}; 
				EndIf
	
				If (t>1 && t<ncol1)
					l1 = newl; Line(l1) = {p1001, p1000};
					l2 = newl; Line(l2) = {p1001, p1002};
					l3 = newl; Line(l3) = {p2001, p2000}; 
					l4 = newl; Line(l4) = {p2001, p2002}; 
		
					l5 = newl; Line(l5) = {p1000, p1003};
					l6 = newl; Line(l6) = {p1003, p1004};
					l7 = newl; Line(l7) = {p1004, p1001}; 
					l8 = newl; Line(l8) = {p2000, p2003}; 
					l9 = newl; Line(l9) = {p2003, p2004}; 
					l10 = newl; Line(l10) = {p2004, p2001}; 
				EndIf
  
				If (t==ncol1)
					l1 = newl; Line(l1) = {p1001, p1000};
					l2 = newl; Line(l2) = {p1001, p25};
					l3 = newl; Line(l3) = {p2001, p2000}; 
					l4 = newl; Line(l4) = {p2001, p53}; 
			
					l5 = newl; Line(l5) = {p1000, p1003};
					l6 = newl; Line(l6) = {p1003, p1004};
					l7 = newl; Line(l7) = {p1004, p1001}; 
					l8 = newl; Line(l8) = {p2000, p2003}; 
					l9 = newl; Line(l9) = {p2003, p2004}; 
					l10 = newl; Line(l10) = {p2004, p2001}; 
				EndIf
				
// 3.- In this section are defined the group lines that will be used later to generate de surfaces of the device  

				thelinesi1[t-1] = l1;
				thelinesi2[t-1] = -l2;
				thelinesd1[t-1] = l3;
				thelinesd2[t-1] = -l4;
	
				thecolumnsi1[t-1] = -l5;
				thecolumnsi2[t-1] = -l6;
				thecolumnsi3[t-1] = -l7;
				thecolumnsi4[t-1] = -l2;
		
				thecolumnsd1[t-1] = l8;
				thecolumnsd2[t-1] = l9;
				thecolumnsd3[t-1] = l10;
				thecolumnsd4[t-1] = l4;
	
				contador = contador+1;
//
			EndFor
//
//
//////////////////////////////////////
//       SURFACES DEFINITION        //
/////////////////////////////////////
//	
// 
// In this section are defined the 3 line loops needed to create the 3 different surfaces of the Microfluidic device.
//
// The Surface 5001 corresponds with the left part of the device, the control channel.
// The Surface 6001 corresponds with the right part of the device, the condition channel.
// The Surface 7001 corresponds with the central part of the device, the gel region.
//
			Line Loop(5000) = {55,2000,200,2002,201,202,-2005,2004,-205,-204,-2003,-203,-2001,thelinesi1[],thelinesi2[],-57, 1001, 103, 1003, 104, 105, -1004, 1005, -102, -101, -1002, -100, -1000};
			Plane Surface(5001) = {5000};
	
			Line Loop(6000) = {3001, 303, 3003, 304, 305, -3004, 3005, -302, -301, -3002, -300, -3000, 56, 4000, 400, 4002, 401, 402, -4005, 4004, -405, -404, -4003, -403, -4001,thelinesd1[],thelinesd2[],-58};
			Plane Surface(6001) = {-6000};

			Line Loop(7000) = {22, 23, 24, 25, 58,thecolumnsd1[],thecolumnsd2[],thecolumnsd3[],thecolumnsd4[],41, 42, 43, 44, -21, -20, -19,thecolumnsi1[],thecolumnsi2[],thecolumnsi3[],thecolumnsi4[],-57, -3, -2, -1};
			Plane Surface(7001) = {-7000};
			
// In this section are defined the colors to assign to the 3 different surfaces crated.
	
			Color Green	{Surface {7001};}
			Color Red {Surface {5001};} 
			Color Blue {Surface {6001};} 
		
		EndIf

	EndIf
//
//////////////////////////////////////
//       PRINT ON COMMAND LINE      //
/////////////////////////////////////
//	
    If (Fabs(p114y)+Fabs(p214y) < (rdep1+rdep2+0.05))
	
		Printf("  Error");
		Printf("  The reservoirs are overlapped for user-defined parameters.");
		
	EndIf	
//	
	If ((hcoltot1+minhnocoltot1) > hdif1)
	
		Printf("  Error");
		Printf("  The posts defined are not compatible with the geometry.");
		
	EndIf
	
