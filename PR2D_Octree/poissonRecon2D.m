function C = poissonRecon2D(ptCloud2d, maxDepth, minDepth, verbose)
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

% TODO: add quadtree
% TODO: qtdecomp, robotics.OccupancyMap3D class, Octree

if nargin < 3
    verbose = false;
end
if maxDepth < minDepth
    disp('maxDepth < minDepth !')
    return;
end

degree = 1;
global dotTable dotdTable ddotdTable
[dotTable, dotdTable, ddotdTable] = getDotTable(degree, minDepth, maxDepth);

time = zeros(5, 1);
tic;
% Create Trees and Samples
minWidth = 2^(-maxDepth);
[pc, T, scale] = normalization(ptCloud2d, 1.5);
pc = pcdownsample2D(pc, minWidth/3);
samples = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);
[tree,samples] = setTree(samples, maxDepth, minDepth);

% Set the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
time(2) = toc() - time(1);
A = setCoefficients(tree);
b = setConstantTerms(tree, samples);
time(3) = toc() - time(2);

% Solve the Linear System
x = A \ b;
time(4) = toc() - time(3);
% TODO: test the influence on speed and effect(ptCloud.Count, depth) of
% scaleFactor, FEM_Basis_dim, weight_Basis_dim, weight_depth, weight_div,
% b_div, grid_div, iso_div, X_div, cgs(). Refer to c++

% Show
if verbose
%     figure
%     quiver(ptCloud2d.Location(:,1), ptCloud2d.Location(:,2), ptCloud2d.Normal(:,1), ptCloud2d.Normal(:,2))
%     title('Input Oriented Points')

%     figure
%     fnplt(bspline(0 : Basis_dim + 1))
%     title('B-Spline')
    
%     figure
%     plot3(samples.Location(:,1), samples.Location(:,2), weight, '.')
%     title('Weight')

%     figure
%     spy(A)
%     title('Coefficients of Linear System')

    figure
    U = tree.width : tree.width : 1;
    [V, U] = meshgrid(U);
    Z = reshape(b, N, N);
    surf(U, V, Z)
    title('Constant Terms of Linear System')
   
    figure
    Z = reshape(x, N, N);
    surf(U, V, Z)
    title('Solution of Linear System')
end

% Extract Contour Line from x
tic;
grid_div = 3;
iso_value = getIsoValue(tree, samples, weight, x, degree, 10);
X = getLinearBsplineSum(tree, x, degree, grid_div, 10);

w = tree.width / grid_div;
U = w : w : 1;
[V, U] = meshgrid(U);
U = (U - 0.5) * scale - T(1);
V = (V - 0.5) * scale - T(2);
Z = reshape(X, grid_div * N, grid_div * N);

% figure, hold on
% plot(ptCloud2d.Location(:,1), ptCloud2d.Location(:,2), '.')
% contour(U, V, Z, [iso_value, iso_value], 'LineWidth', 2);
% hold off

figure
C = contour(U, V, Z, [iso_value, iso_value], 'LineWidth', 2);
title('Isoline')
time(5) = toc();

if verbose
    disp(['Read input into grid:        ',	num2str(time(1))])
    disp(['Got kernel density:          ',	num2str(time(2))])
    disp(['Set FEM constraints:         ',	num2str(time(3))])
    disp(['Linear system solved:        ',	num2str(time(4))])
    disp(['Got piecewise linear curve:  ',	num2str(time(5))])
    disp(' ')
end

end

