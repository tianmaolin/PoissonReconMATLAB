% Shows a couple of sample reconstruction using PR - Poisson Surface
% Reconstruction
%
% Maolin Tian, Tongji University, 2018

ptCloud = pcread('..\data\horse_part.ply');
% ptCloud = pcread('..\data\shpere_uniform_2000.ply');
% ptCloud = pcread('..\data\TjHand1_0600_2100_points_2.ply');
% ptCloud = pcdownsample(ptCloud, 'random',0.5);
minDepth = 4; % max voxel width = 2^-minDepth
maxDepth = 6;
verbose = true;

pcshow(ptCloud)

isoSurface = poissonRecon(ptCloud, minDepth, maxDepth, verbose);





% % make pointCloud of a shpere
% 
% N = 2000;
% 
% % Fibonacci lattice
% n = 1:N;
% phi = (sqrt(5) - 1) / 2;
% z = (2*n - 1) / N - 1;
% x = sqrt(1 - z.*z) .* cos(2*pi*phi*n);
% y = sqrt(1 - z.*z) .* sin(2*pi*phi*n);
% loc = [x; y; z]';
% nor =  - loc;
% shperePtCloud = pointCloud(loc, 'normal', nor);
% % pcshow(shperePtCloud);
% pcwrite(shperePtCloud,'data\shpere_uniform_2000.ply');

% % Spherical coordinate system
% dt = pi / round(sqrt(N));
% a = dt:2*dt:2*pi-dt;
% b = dt:dt:pi-dt;
% [alpha, beta] = meshgrid(a, b);
% x = sin(alpha) .* cos(beta);
% y = sin(alpha) .* sin(beta);
% z = cos(alpha);
% loc = [x(:), y(:), z(:)];
% nor = loc;
% shperePtCloud = pointCloud(loc, 'normal', nor);
% % pcshow(shperePtCloud);
% pcwrite(shperePtCloud,'data\shpere_non_uniform.ply');
