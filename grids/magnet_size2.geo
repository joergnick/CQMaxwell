Point(6) = {-0, 0, -0, 0.1};
Point(7) = {0, 0.1, -0, 0.1};
Point(8) = {0, 0.5, -0, 0.1};
Point(9) = {-0, 0.9, 0, 0.1};
Point(10) = {0, 1, 0, 0.1};
Circle(1) = {7, 8, 9};
Circle(2) = {6, 8, 10};
Line(3) = {6, 7};
Line(4) = {9, 10};
Line Loop(5) = {3, 1, 4, -2};
Ruled Surface(6) = {5};
Translate {0, 0, 2} {
  Duplicata { Surface{6}; }
}
Line(12) = {10, 21};
Line(13) = {9, 17};
Line(14) = {7, 12};
Line(15) = {11, 6};
Line Loop(16) = {12, 11, 15, 2};
Ruled Surface(17) = {16};
Line Loop(18) = {13, -9, -14, 1};
Ruled Surface(19) = {18};
Line Loop(20) = {10, -12, -4, 13};
Plane Surface(21) = {20};
Line Loop(22) = {14, -8, 15, 3};
Plane Surface(23) = {22};
