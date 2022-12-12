
// Points
lc=0.1;
// Base
Point(1) = {0, 0, 0, 0.2*lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
// Top
Point(5) = {0, 0, 1, 0.2*lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {1, 1, 1, lc};
Point(8) = {0, 1, 1, lc};

// Lines
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

// Surfaces
// Base
Line Loop(1) = {-1,-2,-3,-4};
Surface(1) = {1};
// Top
Line Loop(2) = {5:8};
Surface(2) = {2};
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

// Volume
Surface Loop(1) = {1:6};
Volume(1) = {1}; 

// Tagging
Physical Surface("Top")={2};
Physical Surface("Base")={1};
Physical Volume("All")={1};




