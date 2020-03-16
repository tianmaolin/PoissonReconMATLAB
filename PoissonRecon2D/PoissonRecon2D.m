function C = poissonRecon2D(ptCloud2d, minDepth, maxDepth, verbose, curvatureRatio)
%PoissonRecon2D Perform the Poisson Surface Reconstruction algorithm on 2-D
% point cloud. And support quadtree adaptive to curvature.
%
% pointCloud2D object ptCloud2d: Oriented points
% minDepth: max pixel width 2^-minDepth
% maxDepth: min pixel width 2^-maxDepth
% verbose: Display progress information
% curvatureRatio: curvature adaptive ratio to accelarating original Poisson recon
% Contour line object C: boundary of scaned obejct
%
% There is not multigrid method. We use mldivide \ to solve equation here. 
% It uses 3-D Point Cloud Processing introduced in R2015a.
%
% Maolin Tian, 2018

if nargin < 4
    verbose = false;
end
if nargin < 5
    curvatureRatio = 0.1;
end
if maxDepth < minDepth
    error('maxDepth < minDepth!')
end
if maxDepth > 9
    warning('Depth is too big!')
end

% Preprocessing for Dot Product
global valueTable dotTable dotdTable ddotdTable
[valueTable, dotTable, dotdTable, ddotdTable] = valueDotTable(minDepth, maxDepth);

% Create Tree and Samples
time = zeros(5, 1);
tic;
[ptCloud2d, T, scale] = normalization(ptCloud2d, 1.3);
pc = pcdownsample2D(ptCloud2d, 2^(-maxDepth+1));
samp0 = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);
[tree1,samp1] = setTree(samp0, minDepth - 2, maxDepth - 2);

%  Reset (or Refine) Tree if CurvatureRatio is Valid
if curvatureRatio < 1
    weights = getSamplingDensity(samp1, tree1);
    curvatureField = getCurvatureField(samp1, tree1, weights);
    feature = pc.Location(curvatureField < quantile(curvatureField, curvatureRatio),:);
    samp2 = refineTreeAdaptiveToCurvature(ptCloud2d, maxDepth, pc, feature);
    [tree1,samp1] = setTree(samp2, minDepth - 2, maxDepth - 2);
    [tree,samples] = setTree(samp2, minDepth, maxDepth, feature);
else
    [tree,samples] = setTree(samp0, minDepth, maxDepth);
end
time(1) = toc() + time(1);

% Get Sampling Density Weights
tic
weights = getSamplingDensity(samp1, tree1);
time(2) = toc() + time(2);

% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
tic
A = setCoefficients(tree);
b = setConstantTerms(tree, samples, weights);
time(3) = toc();

% Solve the Linear System
% x = cgs(A, b);
x = A \ b;
time(4) = toc() - time(3);

% Show
if verbose
%     figure
%     quiver(ptCloud2d.Location(:,1), ptCloud2d.Location(:,2), ptCloud2d.Normal(:,1), ptCloud2d.Normal(:,2))
%     title('Input Oriented Points'), legend([num2str(ptCloud2d.Count), ' Points'])

%     figure
%     fnplt(bspline(0 : degree + 1))
%     title('B-Spline')
    
%     figure
%     plot3(samples.Location(:,1), samples.Location(:,2), weights, '.')
%     title('Local Density')
%     figure
%     plot3(samples.Location(:,1), samples.Location(:,2), normalWeights, '.')
%     title('Local Average Normal')

    figure, hold on
    plot(tree.center(:,1), tree.center(:,2), '.')
    plot(samples.Location(:,1), samples.Location(:,2), '.')
    legend('tree', 'samples')
    title('Input Points and Tree Center')

    figure
    spy(A)
    title('Coefficients of Linear System')
    legend(['size: ', num2str(size(A,1)), ' * ', num2str(size(A,1))])

    figure, hold on
    truncB = max(abs(quantile(b(b~=0),[0.25, 0.75])));
    plot3(tree.center(:,1), tree.center(:,2), b,'.')
    plot3(tree.center(b<-truncB,1), tree.center(b<-truncB,2), b(b<-truncB),'o')
    plot3(tree.center(b>truncB,1), tree.center(b>truncB,2), b(b>truncB),'*')
    legend('', ['b < ', num2str(-truncB)], ['b > ', num2str(truncB)])
    title('Constant Terms of Linear System')
   
%     figure, hold on
%     res = abs(A * x - b);
%     truncRes1 = quantile(res(res~=0), 0.9);
%     truncRes2 = quantile(res(res~=0), 0.97);
%     plot3(tree.center(:,1), tree.center(:,2), res,'.')
%     plot3(tree.center(res>truncRes1 & res<=truncRes2,1), tree.center(res>truncRes1 & res<=truncRes2,2), res(res>truncRes1 & res<=truncRes2),'*')
%     plot3(tree.center(res>truncRes2,1), tree.center(res>truncRes2,2), res(res>truncRes2),'ro')
%     legend('', ['res > ', num2str(truncRes1)], ['res > ', num2str(truncRes2)])
%     title('Residuals of Linear System')
%     % The residuals of linear system are not significantly correlative with
%     % the reconstruction error, so it is the error of building system that
%     % introduces the actual error. (b in Ax = b is zero at corner because 
%     % the normals there are almost inversely.)
    
%     figure
%     plot3(tree.center(:,1), tree.center(:,2), x,'.')
%     title('Solution of Linear System')
end

% Extract Contour Line from x
tic;
X = basisSum(tree, x);
iso_value = isoValue(tree, samples, x);

w = 2^-maxDepth;
U = w/2:w:1-w/2;
[U,V]= meshgrid(U, U);
Z = griddata(tree.center(:,1), tree.center(:,2), X, U, V, 'linear');
U = double((U - 0.5) * scale - T(1));
V = double((V - 0.5) * scale - T(2));
figure
C = contour(U, V, Z, [iso_value, iso_value], 'LineWidth', 1);
title('Isoline')
time(5) = toc();

if verbose
%     figure, hold on
%     plot3(tree.center(:,1), tree.center(:,2), X,'.')
%     plot3(tree.center(X>iso_value, 1), tree.center(X>iso_value, 2), X(X>iso_value),'*')
%     legend('\chi < isovalue', '\chi > isovalue'), title('\chi')
    
    disp('                             time (s)')
    disp(['Set tree:                    ',  num2str(time(1))])
    disp(['Got weight:                  ',	num2str(time(2))])
    disp(['Set FEM constraints:         ',	num2str(time(3))])
    disp(['Linear system solved:        ',	num2str(time(4))])
    disp(['Got piecewise linear curve:  ',	num2str(time(5))])
%     disp(['Linear system size:        ',	num2str(size(A,1)), ' * ', num2str(size(A,1))])
% Extract isosurface time is o(s^3), because it's not adaptive here.
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ptCloudNormalized, trans, scale] = normalization(ptCloud, scaleFactor)
%normalization Normalize the ptCloud to [0,1]*[0,1]
if nargin < 2
    scaleFactor = 1.25;
end

trans = - [(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2,...
    (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2];
% trans = repmat(trans, ptCloud.Count, 1);
scale = max([ptCloud.XLimits(2) - ptCloud.XLimits(1), ...
    ptCloud.YLimits(2) - ptCloud.YLimits(1)]);
scale = scale * scaleFactor;

location = ptCloud.Location + trans;
location = location / scale + 0.5;
ptCloudNormalized = pointCloud2D(location, ptCloud.Normal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudNormalized = pcdownsample2D(ptCloud, width)
%pcdownsample2D Downsample the 2-D ptCloud

if ptCloud.Count == 0
    ptCloudNormalized = ptCloud;
    return
end

id = ceil(ptCloud.Location / width);
[~ , ia] = unique(id,'rows');

ptCloudNormalized = pointCloud2D(ptCloud.Location(ia,:), ptCloud.Normal(ia,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
