% make pointCloud of a shpere

N = 20000;

% Fibonacci lattice
n = 1:N;
phi = (sqrt(5) - 1) / 2;
z = (2*n - 1) / N - 1;
x = sqrt(1 - z.*z) .* cos(2*pi*phi*n);
y = sqrt(1 - z.*z) .* sin(2*pi*phi*n);
loc = [x; y; z]';
nor = loc;
shperePtCloud = pointCloud(loc, 'normal', nor);
% pcshow(shperePtCloud);
pcwrite(shperePtCloud,'data\shpere_uniform.ply');

% Spherical coordinate system
dt = pi / round(sqrt(N));
a = dt:2*dt:2*pi-dt;
b = dt:dt:pi-dt;
[alpha, beta] = meshgrid(a, b);
x = sin(alpha) .* cos(beta);
y = sin(alpha) .* sin(beta);
z = cos(alpha);
loc = [x(:), y(:), z(:)];
nor = loc;
shperePtCloud = pointCloud(loc, 'normal', nor);
% pcshow(shperePtCloud);
pcwrite(shperePtCloud,'data\shpere_non_uniform.ply');
