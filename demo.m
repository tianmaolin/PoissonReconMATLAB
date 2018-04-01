% Poisson Surface Reconstruction without octree
% Input Oriented points
% Output Triangulated mesh

ptCloud = pcread('data\teapot.ply');
ptCloud = pcdownsample(ptCloud,'random',0.1);
normal = pcnormals(ptCloud);
depth = 5;
% TODO: use griddata() to test
% TODO: Interpolating Gridded Data, Interpolation of 2-D Selections in 3-D
% Grids, Curve Fitting, List of Library Models for Curve and Surface Fitting

% grid = setgrid(ptCloud);
% Normalization
scale = max([ptCloud.XLimits(2) - ptCloud.XLimits(1), ...
    ptCloud.YLimits(2) - ptCloud.YLimits(1), ...
    ptCloud.ZLimits(2) - ptCloud.ZLimits(1)]);
translation = - [(ptCloud.XLimits(2) + ptCloud.XLimits(1)) / 2,...
    (ptCloud.YLimits(2) + ptCloud.YLimits(1)) / 2,...
    (ptCloud.ZLimits(2) + ptCloud.ZLimits(1)) / 2];
translation = repmat(translation, ptCloud.Count, 1);
location = ptCloud.Location;
location = location + translation;
scale = scale * 1.1;
location = location / scale + 0.5;
% change the orientation of teapot
% temp = location(:,1);
% location(:,1) = location(:,3);
% location(:,3) = temp;
pcNormlized = pointCloud(location, 'Normal', normal);

% Discretization
grid.depth = depth;
N = 2^grid.depth;
grid.width = 1 / 2^grid.depth;
grid.size = [1 / grid.width, 1 / grid.width, 1 / grid.width];
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
    c = find(abs(grid.gridijk(:,3) - grid.gridijk(n,3))<= 2);
    grid.neighbor{n} = intersect(a,b);
    grid.neighbor{n} = intersect(grid.neighbor{n},c);
end

% vol = zeros(grid.size);
% for i = 1 : grid.Count
%     vol(grid.gridijk(i,1),grid.gridijk(i,2),grid.gridijk(i,3)) = 1;
% end
% weight = imgaussfilt3(vol, 1);
% for i = 1 : grid.Count
%     grid.weight(i) = weight(grid.gridijk(i,1),grid.gridijk(i,2),grid.gridijk(i,3));
% end


% % [L, v] = SetLaplacianConstraints(grid);
% Vx = zeros(grid.size);
% Vy = zeros(grid.size);
% Vz = zeros(grid.size);
% for i = 1 : grid.Count
%     Vx(grid.gridijk(i,1),grid.gridijk(i,2),grid.gridijk(i,3)) = grid.Normal(i,1) ;
%     Vy(grid.gridijk(i,1),grid.gridijk(i,2),grid.gridijk(i,3)) = grid.Normal(i,2) ;
%     Vz(grid.gridijk(i,1),grid.gridijk(i,2),grid.gridijk(i,3)) = grid.Normal(i,3) ;
% end
% Vx = imgaussfilt3(Vx, 10);
% Vy = imgaussfilt3(Vy, 10);
% Vz = imgaussfilt3(Vz, 10);
% dV = zeros(grid.size);
% dV(2:end,:,:) = (diff(Vx)+diff(Vy)+diff(Vz))/grid.width;
% 
% % \delta X = dV
% n = grid.size(1) + 1;
% I = eye(n,n);
% B = diag(ones(1,n-1),1)+diag(ones(1,n-1),-1)+diag((-2)*ones(1,n));
% A = kron(B,eye(n))+kron(eye(n),B);
% b = dV(:,:,80);
% b = b(:);


% x = L \ v;
% x = A \ b;
% x = reshape(x,n,n);

% iso_value = GetIsoValue(grid);



% GetMCIsoTriangles(grid, iso_value);