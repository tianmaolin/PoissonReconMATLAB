function C = poissonRecon2D(ptCloud2d, minDepth, maxDepth, verbose)
%PoissonRecon2D Perform the Poisson Surface Reconstruction algorithm on 2-D
% point cloud.
%
% pointCloud2D object ptCloud2d: Oriented points
% depth: [2^depth, 2^depth] is grid size
% verbose: Display progress information
% Contour line object C: boundary of scaned obejct
%
% There aren't octree struct and multigrid method. We use pixel grid to
% replace octree, and use mldivide \ to solve the linear system directly.
%
% Maolin Tian, Tongji University, 2018

% TODO: robotics.OccupancyMap3D class, Octree

if nargin < 4
    verbose = false;
end
if maxDepth < minDepth
    disp('maxDepth < minDepth !')
    return;
end

degree = 2;
global valueTable dotTable dotdTable ddotdTable
[valueTable, dotTable, dotdTable, ddotdTable] = getDotTable(degree, minDepth, maxDepth);

time = zeros(5, 1);
tic;
% Create Trees and Samples
[pc, T, scale] = normalization(ptCloud2d, 1.3);
pc = pcdownsample2D(pc, 2^(-maxDepth-1));
samples = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);
[tree,samples] = setTree(samples, minDepth, maxDepth);

time(1) = toc();
% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
weights = getWeight(samples, minDepth - 2 , maxDepth - 2);
time(2) = toc() - time(1);
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

% Show
if verbose
%     figure
%     quiver(ptCloud2d.Location(:,1), ptCloud2d.Location(:,2), ptCloud2d.Normal(:,1), ptCloud2d.Normal(:,2))
%     title('Input Oriented Points'), legend([num2str(ptCloud2d.Count), ' Points'])

%     figure
%     fnplt(bspline(0 : degree + 1))
%     title('B-Spline')
    
%     figure
%     plot3(samples.Location(:,1), samples.Location(:,2), weight, '.')
%     title('Weight')

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
   
%     figure
%     plot3(tree.center(:,1), tree.center(:,2), x,'.')
%     title('Solution of Linear System')
end

tic;
% Extract Contour Line from x
X = basisSum(tree, x);
isoValue = getIsoValue(tree, samples, x);

figure
w = 2^-maxDepth;
U = w/2:w:1-w/2;
[U,V]= meshgrid(U, U);
Z = griddata(tree.center(:,1), tree.center(:,2), X, U, V, 'linear');
U = double((U - 0.5) * scale - T(1));
V = double((V - 0.5) * scale - T(2));
C = contour(U, V, Z, [isoValue, isoValue], 'LineWidth', 1);
title('Isoline')
time(5) = toc();

if verbose
    figure, hold on
    plot3(tree.center(:,1), tree.center(:,2), X,'.')
    plot3(tree.center(X>isoValue, 1), tree.center(X>isoValue, 2), X(X>isoValue),'*')
    legend('\chi < isovalue', '\chi > isovalue'), title('\chi')
    
    disp(['Set tree:        ',          	num2str(time(1))])
    disp(['Got kernel density:          ',	num2str(time(2))])
    disp(['Set FEM constraints:         ',	num2str(time(3))])
    disp(['Linear system solved:        ',	num2str(time(4))])
    disp(['Got piecewise linear curve:  ',	num2str(time(5))])
%     disp(['Linear system size:        ',	num2str(size(A,1)), ' * ', num2str(size(A,1))])
    disp(' ')
   
end

end

