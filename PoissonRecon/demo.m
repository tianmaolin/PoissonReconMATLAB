% Shows a couple of sample reconstruction using PR - Poisson Surface
% Reconstruction
%
% Maolin Tian, Tongji University, 2018

ptCloud = pcread('..\data\shpere_uniform_2000.ply');
% data: shpere_uniform_2000.ply, horse_part.ply, horse.points.ply, bunny.points.ply, eagle.points.ply
% ptCloud = pcread('..\data\TjHand1_0600_2100_points.ply');
ptCloud1 = pcdownsample(ptCloud, 'random', 0.25);
ptCloud2 = pcdownsample(ptCloud, 'random', 0.25);% ptCloud1 intersect ptCloud2 is small enough.
minDepth = 4; % max voxel width = 2^-minDepth
maxDepth = 5;
verbose = true;

% pcshow(ptCloud)

[F,V] = poissonRecon(ptCloud1, minDepth, maxDepth, verbose);

% error
kdOBJ = KDTreeSearcher(V);
[match, mindist] = knnsearch(kdOBJ,ptCloud2.Location);
error = sqrt(mean(mindist.^2));
disp('error =')
disp(error)
% errorPlot = scatter3( ptCloud2.Location(:,1), ptCloud2.Location(:,2), ptCloud2.Location(:,3), [], mindist , '.');
% colormap jet
% color = [mindist/error/2,zeros(length(mindist),1),max(2*error - mindist,0)/error/2];
% cPointCloud = pointCloud(ptCloud2.Location,'Normal',ptCloud2.Normal,'Color',color);
% pcwrite(cPointCloud,'..\data\error.ply');

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
