// Gmsh project created on Tue Jan 14 17:34:38 2020
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, 1, 1, 1.0};
Point(6) = {0, 1, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {1, 5};
Line(3) = {4, 5};
Line(4) = {4, 3};
Delete {
  Line{2};
}
Line(5) = {4, 1};
Line(6) = {3, 2};
Point(7) = {0, 1, 0, 1.0};
Line Loop(7) = {5, 1, -6, -4};
Plane Surface(8) = {7};
Point(8) = {0, 0, 1, 1.0};
Point(9) = {0, 0, 1, 1.0};
Point(10) = {0, 0, 1, 1.0};
Point(11) = {0, 0, 1, 1.0};
Point(12) = {0, 0, 1, 1.0};
Point(13) = {0, 0, 1, 1.0};
Point(14) = {0, 0, 1, 1.0};
Point(15) = {0, 0, 1, 1.0};
Point(16) = {0, 0, 1, 1.0};
Point(17) = {0, 0, 1, 1.0};
Point(18) = {0, 0, 1, 1.0};
Point(19) = {0, 0, 1, 1.0};
Point(20) = {0, 0, 1, 1.0};
Point(21) = {0, 0, 1, 1.0};
Point(22) = {0, 0, 1, 1.0};
Point(23) = {0, 0, 1, 1.0};
Point(24) = {0, 0, 1, 1.0};
Line(9) = {5, 16};
Line(10) = {9, 1};
Point(25) = {1, 1, 1, 1.0};
Point(26) = {1, 0, 1, 1.0};
Line(11) = {25, 26};
Line(12) = {26, 9};
Line(13) = {25, 5};
Delete {
  Surface{8};
}
Plane Surface(14) = {7};
Delete {
  Surface{14};
}
Line Loop(14) = {3, -13, 11, 12, 10, -5};
Plane Surface(15) = {14};
Delete {
  Surface{15};
}
Point(27) = {1, 0, 1, 0.25};
Point(28) = {1, 0, 1, 0.75};
Point(29) = {1, 0.5, 0.7, 0.75};
