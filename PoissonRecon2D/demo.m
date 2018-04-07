ptCloud = ptCloudExample2D(1000);
depth = 5; % grid.size = [2^depth, 2^depth]
Contour = PoissonRecon2D(ptCloud, depth);
x = reshape(Contour, 2^5, 2^5);
surf(x)
