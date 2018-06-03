function C = poissonRecon2D(ptCloud2d, minDepth, maxDepth, verbose)
%PoissonRecon2D Perform the Poisson Surface Reconstruction algorithm on 2-D
% point cloud.
%
% pointCloud2D object ptCloud2d: Oriented points
% minDepth: max pixel width 2^-minDepth
% maxDepth: min pixel width 2^-maxDepth
% verbose: Display progress information
% Contour line object C: boundary of scaned obejct
%
% There is not multigrid method. We use mldivide \ to solve equation. 
% It uses 3-D Point Cloud Processing introduced in R2015a.
%
% Maolin Tian, Tongji University, 2018

% TODO: 效果考虑对比C++版本的
% TODO: 列残差表说明它不是主要误差
% TODO: —Adaptive refinement of the octree based on residuals measured
% at coarser levels, to allow the output mesh complexity to adapt
% not only to sampling density but also to solution quality.


if nargin < 4
    verbose = false;
end
if maxDepth < minDepth
    error('maxDepth < minDepth !')
end

degree = 2;
global valueTable dotTable dotdTable ddotdTable
[valueTable, dotTable, dotdTable, ddotdTable] = valueDotTable(degree, minDepth, maxDepth);

% Create Tree and Samples
time = zeros(5, 1);
tic;
[ptCloud2d, T, scale] = normalization(ptCloud2d, 1.3);
pc = pcdownsample2D(ptCloud2d, 2^(-maxDepth+1));
s0 = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);
[~,s0] = setTree(s0, minDepth, maxDepth);
time(1) = toc();

% Get weights
weights = getLocationWeight(s0, minDepth - 2 , maxDepth - 2);
normalWeights = getNormalWeight(s0, weights, minDepth - 2, maxDepth - 2);
feature = s0.Location(normalWeights < 0.924,:);
% norm([1,0] + [sqrt(2)/2, sqrt(2)/2])/2 = 0.9239 --- 3/4*pi
time(2) = toc() - time(1);

%  Reset ( refine ) tree
tic
location = [];
normal = [];
samW = 2^(-maxDepth+1);
for s = 1:size(feature,1)
    id = ptCloud2d.Location(:,1) < feature(s,1)+samW & ptCloud2d.Location(:,1) > feature(s,1)-samW &...
        ptCloud2d.Location(:,2) < feature(s,2)+samW & ptCloud2d.Location(:,2) > feature(s,2)-samW;
    location = [location; ptCloud2d.Location(id,:)];
    normal = [normal; ptCloud2d.Normal(id,:)];
    id = pc.Location(:,1) < feature(s,1)+samW & pc.Location(:,1) > feature(s,1)-samW &...
        pc.Location(:,2) < feature(s,2)+samW & pc.Location(:,2) > feature(s,2)-samW;
    pc.Location(id,:) = [];
    pc.Normal(id,:) = [];
end
pc2 = pcdownsample2D(pointCloud2D(location,normal), 2^(-maxDepth));
samples = pointCloud2D([pc.Location;pc2.Location], [pc.Normal;pc2.Normal]);
[tree,samples] = setTree(samples, minDepth, maxDepth, feature);
time(1) = toc() + time(1);

% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
tic
weights = getLocationWeight(samples, minDepth - 2 , maxDepth - 2);
time(2) = toc() + time(2);
tic
A = setCoefficients(tree);
b = setConstantTerms(tree, samples, weights);
time(3) = toc();

% Solve the Linear System
% We need refine octree and hanging node to ensure convergence.
% Though I have not found any wrong reconstruction so far without them :).
% x = cgs(A, b);
x = A \ b;
time(4) = toc() - time(3);
% TODO: test the influence on speed and effect(ptCloud.Count, depth) of
% scaleFactor, FEM_Basis_dim, weight_Basis_dim, weight_depth, weight_div,
% b_div, grid_div, iso_div, X_div, \, cgs(). Refer to c++

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

%     figure, hold on
%     plot(tree.center(:,1), tree.center(:,2), '.')
%     plot(samples.Location(:,1), samples.Location(:,2), '.')
%     legend('tree', 'samples')
%     title('Input Points and Tree Center')

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
    
    disp(['Set tree:                    ',  num2str(time(1))])
    disp(['Got weight:                  ',	num2str(time(2))])
    disp(['Set FEM constraints:         ',	num2str(time(3))])
    disp(['Linear system solved:        ',	num2str(time(4))])
    disp(['Got piecewise linear curve:  ',	num2str(time(5)-time(1))])
	% In practice, the time of setting tree does not dominate the actual running time.
%     disp(['Linear system size:        ',	num2str(size(A,1)), ' * ', num2str(size(A,1))])
   
end
    disp(['Total time:                  ',	num2str(sum(time))])
    disp(' ')
    
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

location = [ptCloud.Location,zeros(ptCloud.Count,1)];
normal = [ptCloud.Normal,zeros(ptCloud.Count,1)];

p = pointCloud(location,'Normal',normal);
p = pcdownsample(p,'gridAverage',width);

location = p.Location(:,1:2);
normal = p.Normal(:,1:2);

ptCloudNormalized = pointCloud2D(location, normal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
