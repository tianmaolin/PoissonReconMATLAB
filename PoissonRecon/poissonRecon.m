function [F, V] = poissonRecon(ptCloud, minDepth, maxDepth, ...
  verbose, normalRatio)
%PoissonRecon Perform the Poisson Surface Reconstruction algorithm.
% And support octree adaptive to the normalfield.
%
%    ptCloud:   Oriented points, specified as a pointCloud object.
%   minDepth:   Max grid width = 2^-minDepth
%   maxDepth:   Min grid width = 2^-maxDepth
%    verbose:   show the equation's matrix and reconstruction runtime
%normalRatio:the ratio of high normalfield regions that have max depth
%          F:   the face list of the resulting triangulated mesh
%          V:   the vertex list of the resulting triangulated mesh
%
% There is not multigrid method. We use mldivide \ to solve equation here.
% This function uses 3-D Point Cloud Processing introduced in MATLAB R2015a,
% 3rd part STL\_Export and OcTree (modified) in MathWorks File Exchange.
%
% Maolin Tian, 2018

if nargin < 4
  verbose = false;
elseif nargin < 5
  normalRatio = 0.1;
end
if maxDepth < minDepth
  error('maxDepth < minDepth !')
end
if maxDepth > 7
  warning('Depth is too big!')
end
% TODO: Why can not get a smoothing cube?

% addpath ..\3rdpart\OcTree % OcTreeModified.m
% addpath ..\3rdpart\MarchingCubes
addpath ..\3rdpart\STL_Export

global valueTable dotTable dotdTable ddotdTable
[valueTable, dotTable, dotdTable, ddotdTable] = valueDotTable(minDepth, maxDepth);

% Build octree and samples
% TODO: robotics.OccupancyMap3D class
time = zeros(5, 1);
tic;
[ptCloud, T, scale] = normalization(ptCloud, 1.1);
pc = pcdownsampleConst(ptCloud, 2^(-maxDepth + 1));
samp0 = struct('Count', pc.Count, 'Location', pc.Location, 'Normal', pc.Normal);
[tree1, samp1] = setTree(samp0, minDepth-2, maxDepth-2);
weights = getSamplingDensity(samp1, tree1);

if normalRatio < 1
  normalField = getNormalField(samp1, tree1, weights);
  maxNormW = quantile(normalField, normalRatio);
  feature = samp1.Location(normalField < maxNormW, :);
  [samples, pointsFeature] = refineTreeAdaptiveToNormal(ptCloud, samp1, maxDepth, feature);
  [tree1, samp1] = setTree(samples, minDepth-2, maxDepth-2);
  [tree, samples] = setTree(samples, minDepth, maxDepth, pointsFeature);
else
  [tree, samples] = setTree(samp0, minDepth, maxDepth);
end
time(1) = toc();

% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
weights = getSamplingDensity(samp1, tree1);
time(2) = toc() - time(1);
A = setCoefficients(tree);
b = setConstantTerms(tree, samples, weights);
time(3) = toc() - time(2);

% Solve the Linear System
% x = cgs(A, b);
x = A \ b;
time(4) = toc() - time(3);

% Extract IsoSurface from x
X = basisSum(tree, x);
iso_value = isoValue(tree, samples, x);

w = 2^-maxDepth;
U1 = w / 2:w:1 - w / 2;
[U1, U2, U3] = meshgrid(U1, U1, U1);
Z = griddata(double(tree.center(:, 1)), double(tree.center(:, 2)), double(tree.center(:, 3)) ...
  , X, U1, U2, U3, 'natural');
U1 = double((U1-0.5)*scale-T(1));
U2 = double((U2-0.5)*scale-T(2));
U3 = double((U3-0.5)*scale-T(3));
% % MarchingCubes' quality is not good, though it is fast.
% [F,V] = MarchingCubes(U1, U2, U3, Z, iso_value);
% time(5) = toc();
[F, V] = isosurface(U1, U2, U3, Z, iso_value);
time(5) = toc();

if verbose
  %     figure, hold on
  %     plot3(tree.center(:,1), tree.center(:,2), X,'.')
  %     plot3(tree.center(X>iso_value, 1), tree.center(X>iso_value, 2), X(X>iso_value),'*')
  %     legend('\chi < isovalue', '\chi > isovalue'), title('\chi')
  %     figure, hold on
  %     plot3(samp1.Location(:,1), samp1.Location(:,2), samp1.Location(:,3), '.');
  %     plot3(pc2.Location(:,1), pc2.Location(:,2), pc2.Location(:,3), '.');
  %     title('Feature')

  %     figure
  %     spy(A)
  %     title('Coefficients of Linear System')
  %     legend(['size: ', num2str(size(A,1)), ' * ', num2str(size(A,1))])

  disp(['Set tree:              ', num2str(time(1))])
  disp(['Got kernel density:    ', num2str(time(2))])
  disp(['Set FEM constraints:   ', num2str(time(3))])
  disp(['Linear system solved:  ', num2str(time(4))])
  disp(['Extract isosurface:    ', num2str(time(5))])
end

figure
pa = patch('Faces', F, 'Vertices', V);
isonormals(U1, U2, U3, Z, pa)
pa.FaceColor = 'red';
pa.EdgeColor = 'none';
daspect([1, 1, 1])
view(3);
axis tight
camlight
lighting gouraud
title('IsoSurface')

% solid = ['TotalTime:', num2str(sum(time)-time(5)),...
%     ',MinDepth:', num2str(minDepth), ',MaxDepth:', num2str(maxDepth)];
% STL_Export(V, F, '..\data\recon_result.stl', solid);
disp('Output:')
STL_Export(V, F, char("recon_result.stl"), char("tml10016"));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ptCloudNormalized, trans, scale] = normalization(ptCloud, scaleFactor)
%normalization Normalize the ptCloud to [0,1]*[0,1]
if nargin < 2
  scaleFactor = 1.3;
end

trans = -[(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2, ...
  (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2, ...
  (ptCloud.ZLimits(2) + ptCloud.ZLimits(1)) / 2];
% trans = repmat(trans, ptCloud.Count, 1);
scale = max([ptCloud.XLimits(2) - ptCloud.XLimits(1), ...
  ptCloud.YLimits(2) - ptCloud.YLimits(1), ...
  ptCloud.ZLimits(2) - ptCloud.ZLimits(1)]);
scale = scale * scaleFactor;

location = ptCloud.Location + trans;
location = location / scale + 0.5;
ptCloudNormalized = pointCloud(location, 'Normal', ptCloud.Normal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudOut = pcdownsampleConst(ptCloudIn, width)
%pcdownsampleConst Downsample the 2-D ptCloud, unique points in a box.
% pcdownsample changes data, resulting in error.

if ptCloudIn.Count == 0
  ptCloudOut = ptCloudIn;
  return
end

id = ceil(ptCloudIn.Location/width);
[~, ia] = unique(id, 'rows');

ptCloudOut = pointCloud(ptCloudIn.Location(ia, :), 'Normal', ptCloudIn.Normal(ia, :));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
