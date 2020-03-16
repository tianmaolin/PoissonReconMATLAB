% Show some surface reconstruction using Poisson Surface Reconstruction
%
% Maolin Tian, 2018

% Data: shpere_non_uniform.points, horse.points, bunny.points,
% Armadillo.points, dragon.points, eagle.points.
ptCloud = pcread('..\data\horse.points.ply');
minDepth = 4; % max voxel width = 2^-minDepth
maxDepth = 6; % maxDepth should < 8
verbose = false;

ptCloudTrain = pointCloud(ptCloud.Location(1:2:end,:), 'Normal', ptCloud.Normal(1:2:end,:));
ptCloudTest = pointCloud(ptCloud.Location(2:2:end,:), 'Normal', ptCloud.Normal(2:2:end,:));
pcshow(ptCloudTrain)
title('Input Point Cloud')

% Poisson Surface Reconstruction
[F,V] = poissonRecon(ptCloudTrain, minDepth, maxDepth, verbose);

% % Error
% [error, dist] = getError(V, ptCloudTest.Location);
% % figure, axis equal, hold on
% % % scatter3( V(:,1), V(:,2), V(:,3),'.');
% % scatter3( ptCloudTest.Location(dist>2*error,1), ptCloudTest.Location(dist>2*error,2), ptCloudTest.Location(dist>2*error,3),'.');
% % errorPlot = scatter3( ptCloudTest.Location(:,1), ptCloudTest.Location(:,2), ptCloudTest.Location(:,3), [], dist , '.');
% % colormap jet
% % color = [dist/error/2,zeros(length(dist),1),max(2*error - dist,0)/error/2];
% % cPointCloud = pointCloud(ptCloudTest.Location,'Color',color, 'Normal', ptCloudTest.Normal);
% % pcwrite(cPointCloud,'..\data\error.ply');
