// Domain Constants
// Core
Cxmin = -1000.;    // min x
Cxmax =  1000.;    // max x
Cymin = -1500.;    // min y
Cymax =  1500.;    // max y
Czmin = -1500.;    // min z
GrndL =  0.;       // Ground Level

// Buffer
Bxmin = Cxmin - 1000.;   // min x
Bxmax = Cxmax + 1000.;   // max x
Bymin = Cymin - 1000.;   // min y
Bymax = Cymax + 1000.;   // max y
Bzmin = Czmin - 1000.;   // min z

// Air Layer
airHt = 1000.;

// Mesh element sizes
//Core
MCground = 100.;     // at ground level
MCbase = 200.;       // at base of core

// Buffer
MBbase = 600.;
MBground = 300.;

// Air
MAir = 600.;

// sensor mesh length
MSens = 10;

//------------------------------------------------------
// Points

// CORE
// Base
Point(1)  = {Cxmin, Cymin, Czmin, MCbase};
Point(2)  = {Cxmax, Cymin, Czmin, MCbase};
Point(3)  = {Cxmax, Cymax, Czmin, MCbase};
Point(4)  = {Cxmin, Cymax, Czmin, MCbase};
// Ground
Point(5)  = {Cxmin, Cymin, 0., MCground};
Point(6)  = {Cxmax, Cymin, 0., MCground};
Point(7)  = {Cxmax, Cymax, 0., MCground};
Point(8)  = {Cxmin, Cymax, 0., MCground};

// BUFFER
// Base
Point(9)  = {Bxmin, Bymin, Bzmin, MBbase};
Point(10)  = {Bxmax, Bymin, Bzmin, MBbase};
Point(11)  = {Bxmax, Bymax, Bzmin, MBbase};
Point(12)  = {Bxmin, Bymax, Bzmin, MBbase};
// Ground
Point(13)  = {Bxmin, Bymin, 0., MBground};
Point(14)  = {Bxmax, Bymin, 0., MBground};
Point(15)  = {Bxmax, Bymax, 0., MBground};
Point(16)  = {Bxmin, Bymax, 0., MBground};

// AIR Layer
Point(17)  = {Bxmin, Bymin, airHt, MAir};
Point(18)  = {Bxmax, Bymin, airHt, MAir};
Point(19)  = {Bxmax, Bymax, airHt, MAir};
Point(20)  = {Bxmin, Bymax, airHt, MAir};


//------------------------------------------------------
// Lines

// CORE
// Base
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Top
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
// Sides
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// BUFFER
Line(13) = {9, 10};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};
// Top
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};
// Sides
Line(21) = {9, 13};
Line(22) = {10, 14};
Line(23) = {11, 15};
Line(24) = {12, 16};

// AIR
// Top
Line(25) = {17, 18};
Line(26) = {18, 19};
Line(27) = {19, 20};
Line(28) = {20, 17};
// Sides
Line(29) = {13, 17};
Line(30) = {14, 18};
Line(31) = {15, 19};
Line(32) = {16, 20};

//------------------------------------------------------
// Surfaces
// all horizontal surfaces face down for z < 0 and up for z>=0
// CORE
// core horizontal surfaces
// Surfaces
// Base
Line Loop(1) = {-1,-2,-3,-4};
Surface(1) = {1};
Physical Surface("coreBase") = {1};
// Top
Line Loop(2) = {5:8};
Surface(2) = {2};
Physical Surface("coreInterface") = {2};
// Front
Line Loop(3) = {1,10,-5,-9};
Surface(3) = {3};
// Back
Line Loop(4) = {3, 12, -7, -11};
Surface(4) = {4};
// Left
Line Loop(5) = {4, 9, -8, -12};
Surface(5) = {5};
// Right
Line Loop(6) = {2, 11,-6, -10};
Surface(6) = {6};

// BUFFER
// Base
Line Loop(7) = {-13,-14,-15,-16};
Surface(7) = {7};
Physical Surface("BufferBase") = {7};
// Top
Line Loop(8) = {17:20};
Surface(8) = {8,-2};
Physical Surface("BufferInterface") = {8};
// Front
Line Loop(9) = {13, 22,-17,-21};
Surface(9) = {9};
// Back
Line Loop(10) = {15, 24, -19, -23};
Surface(10) = {10};
// Left
Line Loop(11) = {16, 21, -20, -24};
Surface(11) = {11};
// Right
Line Loop(12) = {14, 23, -18, -22};
Surface(12) = {12};

// AIR
// Top
Line Loop(13) = {25:28};
Surface(13) = {13};
Physical Surface("Top") = {13};
// Front
Line Loop(14) = {17, 30, -25, -29};
Surface(14) = {14};
// Back
Line Loop(15) = {19, 32, -27, -31};
Surface(15) = {15};
// Left
Line Loop(16) = {20, 29, -28, -32};
Surface(16) = {16};
// Right
Line Loop(17) = {18, 31, -26, -30};
Surface(17) = {17};

//  SENSOR LOCATIONS
k=newp;
Point(k)={10.,-100., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={100., 100., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={150., -500., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={500.,-100., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={-500., 120, 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={150.,-600., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={300., -700., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={100., 800., 0.,MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={150., 650., 0., MSens};
Point{k} In Surface{2};
k=newp;
Point(k)={-300., 1000., 0., MSens};
Point{k} In Surface{2};

// VOLUMES
// core 
Surface Loop(18) = {1:6};
Volume(1) = {18}; 
Physical Volume("coreGround") = {1};
// buffer 
Surface Loop(19) = {7:12,-1,-3,-4,-5,-6};
Volume(2) = {19};
Physical Volume("buffer") = {2};
// air
Surface Loop(20) = {2,8,13:17};
Volume(3) = {20}; 
Physical Volume("bufferAir") = {3};
