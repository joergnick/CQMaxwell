// Gmsh project created on Tue Sep 22 15:23:11 2020
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Circle(1) = {2, 1, 2};
Point(3) = {1.1, 0, 0, 1.0};
Point(4) = {1.1, 0, 0.1, 1.0};
Point(5) = {1.2, 0, 0, 1.0};
Point(6) = {1.1, 0, -0.1, 1.0};
Circle(2) = {2, 3, 6};
Circle(3)= { 6,3,5};
Circle(4)={5,3,4};
Circle(5)={4,3,2};

