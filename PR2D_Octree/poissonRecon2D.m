function C = poissonRecon2D(ptCloud2d, depth, verbose)
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

time = zeros(5, 1);
tic;
% Create Grids and Samples
N = 2^depth;
grid = struct('depth', depth, 'width', 1/N, 'size', [N, N]);
[pc, T, scale] = normalization(ptCloud2d, 1.3);
samples = struct('Count', pc.Count, 'Location', pc.Location,'Normal', pc.Normal);

% Create Maps between Grid and Samples
% grid_ind is the top-right corner of pixel grid containing the sample
samples.grid_ind = ceil(samples.Location / grid.width);
% convert (i,j) to n = i + (j - 1) * N
samples.grid_ind = samples.grid_ind(:, 1) + (samples.grid_ind(:, 2) - 1) * N;
% sample_ind is the cell of samples' indices contained by the pixel grid
grid.sample_ind = cell(N * N, 1);
for n = 1 : N * N
    grid.sample_ind{n} = find(samples.grid_ind == n);
end
time(1) = toc();

% Get the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
Basis_dim = 2;
weight = getWeight(grid, samples, grid.depth - 1, Basis_dim, 20);
time(2) = toc() - time(1);
A = getCoefficients(grid, Basis_dim);
b = getConstantTerms(grid, samples, weight, Basis_dim, 5);
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

    figure
    fnplt(bspline(0 : Basis_dim + 1))
    title('B-Spline')
    
%     figure
%     plot3(samples.Location(:,1), samples.Location(:,2), weight, '.')
%     title('Weight')

    figure
    spy(A)
    title('Coefficients of Linear System')

    figure
    U = grid.width : grid.width : 1;
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
iso_value = getIsoValue(grid, samples, weight, x, Basis_dim, 10);
X = getLinearBsplineSum(grid, x, Basis_dim, grid_div, 10);

w = grid.width / grid_div;
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
