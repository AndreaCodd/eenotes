Cxmin =-2.000000e+04;
Cxmax =2.000000e+04;
Cymin =-5.000000e+04;
Cymax =0.000000e+00;
cg =1.000000e+02;
cb =1.000000e+03;
Bx = 1.000000e+05;
By = 1.000000e+05;
Bair = 1.000000e+05;
bb = 5.000000e+04;
bg = 1.000000e+04;
ba = 5.000000e+04;
snum = 1.500000e+01;
sspan = 1.800000e+04;
sg = 6.428571e+01;

sx= sspan/(snum-1);
s1= sspan/2.;


// 2D
// CORE
// core Points

Point(1) = {Cxmin, Cymin, 0, cb};
Point(2) = {Cxmax, Cymin, 0, cb};
Point(3) = {Cxmax, Cymax, 0, cg};
Point(4) = {Cxmin, Cymax, 0, cg};

Point(5) = {Cxmin - Bx , Cymin-By, 0.0, bb};
Point(6) = {Cxmax + Bx , Cymin-By, 0.0, bb};
Point(7) = {Cxmax + Bx , 0.0, 0.0, bg};
Point(8) = {Cxmax + Bx , Bair, 0.0, ba};
Point(9) = {Cxmin - Bx , Bair, 0.0, ba};
Point(10) = {Cxmin - Bx , 0.0, 0.0, bg};

// mesh nodes 
k=newp;
For i In {0:snum-1}
    Point(k+i)={ s1, 0.0, 0.0, sg};
    s1=s1-sx;
EndFor

// core lines, ground, bottom, sides
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,k};
For i In {1:snum-1}
    Line(3+i) = {k+i-1,k+i};
EndFor
l=newl;
Line(l) = {k+snum-1,4};
Line(l+1) = {4,1};

Line Loop(1) = {1:l+1};
Plane Surface(1) = {1};
Physical Surface("core") = {1};

m=newl;
Line(m) = {5,6};
Line(m+1) = {6,7};
Line(m+2) = {7,8};
Line(m+3) = {8,9};
Line(m+4) = {9,10};
Line(m+5) = {10,5};
Line(m+6) = {10,4};
Line(m+7) = {3,7};

Line Loop(2) = {m,m+1,-m-7,-2,-1,-m+1,-m-6,m+5};
Plane Surface(2) = {2};
Physical Surface("buffer") = {2};

Line Loop(3) = {-m-7,3:l,-m-6,-m-4,-m-3,-m-2};
Plane Surface(3) = {3};
Physical Surface("air") = {3};
