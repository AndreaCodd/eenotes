
// Points
lc=0.1;
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, 0.2*lc};
Point(4) = {0, 1, 0, 0.1*lc};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Surface
Line Loop(1) = {1:4};
Surface(1) = {1};

// Tagging
Physical Line("Top")={3};
Physical Line("Base")={1};
Physical Surface("All")={1};
