function [F,V] = poissonRecon(ptCloud, minDepth, maxDepth, verbose)
%PoissonRecon Perform the Poisson Surface Reconstruction algorithm.
%
% pointCloud:   Oriented points, specified as a pointCloud object.
%   minDepth:   Max grid width = 2^-minDepth
%   maxDepth:   Min grid width = 2^-maxDepth
%    verbose:   show the equation's matrix and reconstruction runtime
%          F:   the face list of the resulting triangulated mesh
%          V:   the vertex list of the resulting triangulated mesh
%
% There is not multigrid method. We use mldivide \ to solve equation. 
% It uses 3-D Point Cloud Processing introduced in MATLAB R2015a.
%
% Maolin Tian, Tongji University, 2018

if nargin < 4
    verbose = false;
end
if maxDepth < minDepth
    error('maxDepth < minDepth !')
end

% addpath ..\3rdpart\OcTree % OcTreeModified.m
% addpath ..\3rdpart\MarchingCubes
addpath ..\3rdpart\STL_Export

degree = 2;
global valueTable dotTable dotdTable ddotdTable
[valueTable, dotTable, dotdTable, ddotdTable] = valueDotTable(degree, minDepth, maxDepth);

% Create octree and samples
% TODO: robotics.OccupancyMap3D class
time = zeros(5, 1);
tic;
[ptCloud, T, scale] = normalization(ptCloud, 1.2);
pc = pcdownsample(ptCloud,'gridAverage', 2^(-maxDepth+1));
samp0 = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);
[tree0,samp0] = setTree(samp0, minDepth - 2, maxDepth - 2);
time(1) = toc();

% Get weights
maxNormW = 0.708;
weights = getLocationWeight(samp0, tree0);
normalWeights = getNormalWeight(samp0, tree0, weights);
feature = samp0.Location(normalWeights < maxNormW,:);
% norm([1,0] + [sqrt(2)/2, sqrt(2)/2])/2 = 0.9239 --- 3/4*pi
time(2) = toc() - time(1);

%  Reset ( refine ) tree
tic
location = [];
normal = [];
samW = 2^-maxDepth;
for s = 1:size(feature,1)
    id = ptCloud.Location(:,1) < feature(s,1)+samW & ptCloud.Location(:,1) > feature(s,1)-samW &...
        ptCloud.Location(:,2) < feature(s,2)+samW & ptCloud.Location(:,2) > feature(s,2)-samW &...
        ptCloud.Location(:,3) < feature(s,3)+samW & ptCloud.Location(:,3) > feature(s,3)-samW;
    location = [location; ptCloud.Location(id,:)];
    normal = [normal; ptCloud.Normal(id,:)];
    id = samp0.Location(:,1) < feature(s,1)+samW & samp0.Location(:,1) > feature(s,1)-samW &...
        samp0.Location(:,2) < feature(s,2)+samW & samp0.Location(:,2) > feature(s,2)-samW &...
        samp0.Location(:,3) < feature(s,3)+samW & samp0.Location(:,3) > feature(s,3)-samW;
    samp0.Location(id,:) = [];
    samp0.Normal(id,:) = [];
end
pc2 = pcdownsample(pointCloud(location, 'Normal', normal), 'gridAverage', 2^(-maxDepth));
samples = struct('Count', size(samp0.Location,1) + size(pc2.Location,1), 'Location', [samp0.Location;pc2.Location],'Normal', [samp0.Normal;pc2.Normal]);
pointsFeature = [false(size(samp0.Location,1),1); true(size(pc2.Location,1),1)];
[tree,samples] = setTree(samples, minDepth, maxDepth, pointsFeature);
[tree1,samp1] = setTree(samples, minDepth - 2, maxDepth - 2);
time(1) = toc() + time(1);

% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
tic
weights = getLocationWeight(samp1, tree1);
time(2) = toc() + time(2);
tic
A = setCoefficients(tree);
b = setConstantTerms(tree, samples, weights);
time(3) = toc() - time(2);

% Solve the Linear System
% x = cgs(A, b);
x = A \ b;
time(4) = toc() - time(3);
% TODO: test the influence on speed and effect(ptCloud.Count, depth) of
% scaleFactor, FEM_Basis_dim, weight_Basis_dim, weight_depth, weight_div,
% b_div, grid_div, iso_div, X_div, \, cgs(). Refer to c++

% Extract IsoSurface from x
X = basisSum(tree, x);
iso_value = isoValue(tree, samples, x);

w = 2^-maxDepth;
U1 = w/2:w:1-w/2;
[U1,U2,U3]= meshgrid(U1, U1, U1);
Z = griddata(double(tree.center(:,1)), double(tree.center(:,2)), double(tree.center(:,3))...
    , X, U1, U2, U3, 'linear');
U1 = double((U1 - 0.5) * scale - T(1));
U2 = double((U2 - 0.5) * scale - T(2));
U3 = double((U3 - 0.5) * scale - T(3));
% % MarchingCubes' quality is not good, though it is fast.
% [F,V] = MarchingCubes(U1, U2, U3, Z, iso_value);
% time(5) = toc() - time(4);
[F, V] = isosurface(U1, U2, U3, Z, iso_value, 'verbose');
time(5) = toc();

if verbose
%     figure, hold on
%     plot3(tree.center(:,1), tree.center(:,2), X,'.')
%     plot3(tree.center(X>iso_value, 1), tree.center(X>iso_value, 2), X(X>iso_value),'*')
%     legend('\chi < isovalue', '\chi > isovalue'), title('\chi')
    
    figure
    spy(A)
    title('Coefficients of Linear System')
    legend(['size: ', num2str(size(A,1)), ' * ', num2str(size(A,1))])

    figure
    p = patch('Faces',F,'Vertices',V);
    isonormals(U1, U2, U3, Z, p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3);
    axis tight
    camlight
    lighting gouraud
    title('IsoSurface')

    disp(['Set tree:              ',    num2str(time(1))])
    disp(['Got kernel density:    ',	num2str(time(2))])
    disp(['Set FEM constraints:   ',	num2str(time(3))])
    disp(['Linear system solved:  ',	num2str(time(4))])
    disp(['Extract isosurface:    ',	num2str(time(5))])
%     disp(['Linear system size:        ',	num2str(size(A,1)), ' * ', num2str(size(A,1))])
   
end
disp(['Total Time:    ',	num2str(sum(time)-time(1))])
disp(' ')
solid = ['TotalTime:', num2str(sum(time)-time(1)),...
    ',MinDepth:', num2str(minDepth), ',MaxDepth:', num2str(maxDepth)];
STL_Export(V, F, '..\data\recon_result.stl', solid);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ptCloudNormalized, trans, scale] = normalization(ptCloud, scaleFactor)
%normalization Normalize the ptCloud to [0,1]*[0,1]
if nargin < 2
    scaleFactor = 1.3;
end

trans = - [(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2,...
    (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2,...
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
