% Shows a couple of sample reconstruction using PR - Poisson Surface
% Reconstruction
%
% Maolin Tian, Tongji University, 2018

ptCloud = pcread('..\data\horse.points.ply');
% data: shpere_uniform.ply, horse.points.ply, bunny.points.ply, eagle.points.ply
% ptCloud = pcread('..\data\TjHand1_0600_2100_points.ply');
ptCloud1 = pointCloud(ptCloud.Location(1:2:end,:), 'Normal', ptCloud.Normal(1:2:end,:));
ptCloud2 = pointCloud(ptCloud.Location(2:2:end,:), 'Normal', ptCloud.Normal(2:2:end,:));
minDepth = 4; % max voxel width = 2^-minDepth
maxDepth = 5;
verbose = true;

% pcshow(ptCloud)

[F,V] = poissonRecon(ptCloud1, minDepth, maxDepth, verbose);

% error
[error, dist] = getError(V, ptCloud2.Location);
figure, axis equal, hold on
% scatter3( V(:,1), V(:,2), V(:,3),'.');
scatter3( ptCloud2.Location(dist>1.5*error,1), ptCloud2.Location(dist>1.5*error,2), ptCloud2.Location(dist>1.5*error,3),'.');
% errorPlot = scatter3( V(:,1), V(:,2), V(:,3), [], dist , '.');
% colormap jet
% color = [dist/error/2,zeros(length(dist),1),max(2*error - dist,0)/error/2];
% cPointCloud = pointCloud(P2,'Color',color);
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
