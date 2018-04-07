function C = PoissonRecon2D(ptCloud2d, depth)
%PoissonRecon2D 2-D Poisson Surface Reconstruction
% pointCloud2D object ptCloud2d: Oriented points
% depth: [2^depth, 2^depth] is grid size
% Contour object C: boundary of scaned obejct
% 
% There aren't octree struct and multigrid method. We use pixel grid to
% replace octree, and use backslash operator \ to solve the linear system directly.
% Maolin Tian, Tongji University, 2018

% Create Grids and Samples
N = 2^depth;
grid = struct('depth', depth, 'width', 1/N, 'size', [N, N]);
pc = normalization(ptCloud2d, 1.5);
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

% Get the FEM Coefficients and Constant Terms
% Paper: Kazhdan, Bolitho, and Hoppe. Poisson Surface Reconstruction. 2006
A = getCoefficients(grid);
weight = getWeight(grid, samples, grid.depth - 1, 100);
b = getConstantTerms(grid, samples, weight, 5);

% Solve the Linear System
x = A \ b;
% x = cgs(A,b)

% Show the Results
% TODO: quiver(pc),spy(A),plot3(weight),surf(b),surf(X)
% TODO: show basic function
% TODO: iso_value = GetIsoValue(grid);
% TODO: GetMCIsoTriangles(grid, iso_value);
% X = grid.width:grid.width:1;
% [Y, X] = meshgrid(X);
% iso_value = 1e-5:5e-6:2e-5;
% C = contour(X,Y,x,iso_value);
C = x;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudNormalized = normalization(ptCloud, scaleFactor)
%normalization Normalize the ptCloud to [0,1]*[0,1]
if nargin < 2
    scaleFactor = 1.25;
end

T = - [(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2,...
                 (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2];
T = repmat(T, ptCloud.Count, 1);
scale = max([ptCloud.XLimits(2) - ptCloud.XLimits(1), ...
             ptCloud.YLimits(2) - ptCloud.YLimits(1)]);
scale = scale * scaleFactor;

location = ptCloud.Location + T;
location = location / scale + 0.5;
ptCloudNormalized = pointCloud2D(location, ptCloud.Normal);
end
