// Box Corner
x_box = DefineNumber[ 0, Name "Parameters/box/x_box" ];
y_box = DefineNumber[ 0, Name "Parameters/box/y_box" ];
z_box = DefineNumber[ 0, Name "Parameters/box/z_box" ];

// Box Dimensions
L = DefineNumber[ 1, Name "Parameters/box/L" ];
W = DefineNumber[ 1, Name "Parameters/box/W" ];
H = DefineNumber[ 1, Name "Parameters/box/H" ];

// Refinement Parameters
ref_low = DefineNumber[ 0.1, Name "Parameters/refinement/low" ];
ref_high = DefineNumber[ 0.01, Name "Parameters/refinement/high" ];

// Line of High Refinement
x0 = DefineNumber[ x_box + 0.7*L, Name "Parameters/line/x0" ];
y0 = DefineNumber[ y_box + 0.7*W, Name "Parameters/line/y0" ];
z0 = DefineNumber[ z_box, Name "Parameters/line/z0" ];

x1 = DefineNumber[ x_box, Name "Parameters/line/x1" ];
y1 = DefineNumber[ y_box + 0.3*W, Name "Parameters/line/y1" ];
z1 = DefineNumber[ z_box + 0.3*H, Name "Parameters/line/z1" ];

Point(9) =  {x0, y0, z0, ref_high};
Point(10) = {x1, y1, z1, ref_high};
Line(13) = {9, 10};

// Define Box Points
Point(1) = {x_box + 0, y_box + 0, z_box + 0, ref_low};
Point(2) = {x_box + 0, y_box + W, z_box + 0, ref_low};
Point(3) = {x_box + L, y_box + W, z_box + 0, ref_low};
Point(4) = {x_box + L, y_box + 0, z_box + 0, ref_low};
Point(5) = {x_box + 0, y_box + 0, z_box + H, ref_low};
Point(6) = {x_box + 0, y_box + W, z_box + H, ref_low};
Point(7) = {x_box + L, y_box + W, z_box + H, ref_low};
Point(8) = {x_box + L, y_box + 0, z_box + H, ref_low};

// Define Box Lines
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 5};
Line(8) = {5, 6};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Define Box Loops and Surfaces
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {9, -7, -12, 3};
Curve Loop(4) = {9, 8, -10, -4};
Curve Loop(5) = {10, 5, -11, -1};
Curve Loop(6) = {11, 6, -12, -2};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

// Define Surface Loop and Volume
Surface Loop(1) = {1, -2, -3, 4, 5, 6};
Volume(1) = {1};

// Place Line of Refinement in Volume
Line{13} In Volume{1};
Point{9} In Surface{1};
Point{10} In Surface{4};
