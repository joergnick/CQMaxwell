cl__1 = 1;
Point(1) = {1, 0, 0, 1};
Point(2) = {0.9, 0, -0, 1};
Point(3) = {0.5, 0.1, -0, 1};
Point(4) = {0.5, -0, -0, 1};
Point(5) = {0.1, 0, -0, 1};
Point(6) = {-0, 0, -0, 1};
Line(1) = {6, 5};
Line(2) = {1, 2};
Circle(3) = {5, 4, 2};
Delete {
  Line{3};
}
Circle(3) = {2, 4, 5};
Circle(4) = {1, 4, 6};
Line Loop(5) = {1, -3, -2, 4};
Ruled Surface(6) = {5};
Translate {0, 0, 1} {
  Duplicata { Surface{6}; }
}


Rotate {{0, 0, 0}, {0, 0, 1}, Pi/4} {
  Duplicata { Surface{7}; }
}
Rotate {{0, 0, 0}, {0, 0, 1}, Pi/4} {
  Duplicata { Surface{7}; }
}
Line(15) = {7, 6};
Line(16) = {1, 17};
Delete {
  Line{16};
}
Line(16) = {17, 1};
Line Loop(17) = {15, -4, -16, 11};
