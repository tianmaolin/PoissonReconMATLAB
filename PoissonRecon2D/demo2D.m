% 2D Poisson Surface Reconstruction without octree
% Input Oriented points
% Output Piecewise linear boundary

ptCloud = circle2D(1000);
depth = 6;

% Normalization
scale = max([ptCloud.XLimits(2) - ptCloud.XLimits(1), ...
    ptCloud.YLimits(2) - ptCloud.YLimits(1)]);
translation = - [(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2,...
    (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2];
translation = repmat(translation, ptCloud.Count, 1);
location = ptCloud.Location;
location = location + translation;
scale = scale * 1.25;
location = location / scale + 0.5;
pcNormlized = pointCloud2D(location, ptCloud.Normal);

% Discretization
grid.depth = depth;
N = 2^grid.depth;
grid.width = 1 / 2^grid.depth;
grid.size = [1 / grid.width, 1 / grid.width];
grid.Count = pcNormlized.Count;
grid.Location = pcNormlized.Location;
grid.Normal = pcNormlized.Normal;
% gridijk is the top-right corner of grid
grid.gridijk = ceil(pcNormlized.Location / grid.width);
grid.gridLocation = grid.gridijk * grid.width;
% quadrant() = 1 means that point is in top/right, else = 0
grid.quadrant = (grid.gridLocation - grid.Location < grid.width/2);
grid.neighbor = cell(grid.Count,1);
for n = 1:grid.Count
    a = find(abs(grid.gridijk(:,1) - grid.gridijk(n,1))<= 2);
    b = find(abs(grid.gridijk(:,2) - grid.gridijk(n,2))<= 2);
    grid.neighbor{n} = intersect(a,b);
end

% [A, b] = SetLaplacianConstraints(grid);
A = getSystem(grid.width);
weight = getWeight(grid, grid.depth - 1, 0.01);
b = getb(grid, weight, 5);

% solve x = A \ b
x = A \ b;

% iso_value = GetIsoValue(grid);

% GetMCIsoTriangles(grid, iso_value);