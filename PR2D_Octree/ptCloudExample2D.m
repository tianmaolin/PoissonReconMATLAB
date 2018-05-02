function ptCloud2D = ptCloudExample2D(model, N)
%ptCloud2D make 2D pointCloud examples
%
% Maolin Tian, Tongji University, 2018

if model == 'circle'
dt = pi / round(N);
% alpha = [0:2*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:dt:2*pi-dt];
alpha = 0:2*dt:2*pi-dt;
x = cos(alpha);
y = sin(alpha);
location = [x(:), y(:)];
normal = - [x(:), y(:)];
ptCloud2D = pointCloud2D(location, normal);
return
end

if strcmp(model,'triangle')
    x0 = 0:1/N:1;
    y0 = 1 - x0;
    x = [x0, x0, zeros(1, length(x0))];
    y = [zeros(1, length(x0)), y0, y0];
    location = [x(:), y(:)];
    n1 = repmat([0,1], length(x0),1);
    n2 = repmat([-sqrt(2)/2, -sqrt(2)/2],length(x0),1);
    n3 = repmat([1,0], length(x0),1);
    normal = [n1; n2; n3];
    ptCloud2D = pointCloud2D(location, normal);
    ptCloud2D = AddNoise(ptCloud2D, 0.001);
    
elseif strcmp(model,'circle')
    dt = pi / round(N);
%     alpha = [0:3*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:10*dt:2*pi-dt];
    alpha = [0:dt:pi/2-dt, pi/2:4*dt:3/2*pi-dt, 3/2*pi:10*dt:13/8*pi];
    % alpha = 0:2*dt:2*pi-dt;
    x = cos(alpha);
    y = sin(alpha);
    location = [x(:), y(:)];
    normal = - [x(:), y(:)];
    ptCloud2D = pointCloud2D(location, normal);
    ptCloud2D = AddNoise(ptCloud2D, 0.05);
    
elseif strcmp(model,'Armadillo')
    % Author: The Stanford 3D Scanning Repository
    ptCloud3D = pcread('..\data\Armadillo.ply');
    XZ_i = find(abs(ptCloud3D.Location(:,1)) < 1);
    location = ptCloud3D.Location(XZ_i,:);
    location = [- location(:, 3), location(:, 2)];
    normal = ptCloud3D.Normal(XZ_i,:);
    normal = [- normal(:, 3), normal(:, 2)];
    normal = normal./(normal(:,1).^2 + normal(:,2).^2);
    ptCloud2D = pointCloud2D(location, normal);
%     ptCloud2D = AddNoise(ptCloud2D, 0.01);
   
elseif strcmp(model,'Dragon')
    % Author: The Stanford 3D Scanning Repository
    ptCloud3D = pcread('..\data\dragon.ply');
    XZ_i = find(abs(ptCloud3D.Location(:,3) + 0.0045) < 0.0005 ...
        & ptCloud3D.Location(:,1) < -0.015 ...
        & ptCloud3D.Location(:,1) > -0.105 ...
        & ptCloud3D.Location(:,2) < 0.19 ...
        & ptCloud3D.Location(:,2) > 0.12);
    location = ptCloud3D.Location(XZ_i,:);
    location = [location(:, 1), location(:, 2)];
    normal = ptCloud3D.Normal(XZ_i,:);
    normal = [normal(:, 1), normal(:, 2)];
    normal = normal./(normal(:,1).^2 + normal(:,2).^2);
    ptCloud2D = pointCloud2D(location, normal);
%     ptCloud2D = pcdownsample2D(ptCloud2D, 0.05);
   
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pcNoisy = AddNoise(ptCloud2D, e)
%AddNoise add noise to ptCloud2D
% loc_e = loc + e * e1 * n + e/2 * e2 * t, normal_e = normal + e/2 * e3
% e1, e2, e3 ~ Norm(0,1). t is vertical to normal.
N = ptCloud2D.Count;
e1 = normrnd(0,  e, [N,1]);
e2 = normrnd(0, e/2, [N,1]);
e3 = normrnd(0, e/2, [N,1]);
ptCloud2D.Location = ptCloud2D.Location + e1 .* ptCloud2D.Normal;
t = [-ptCloud2D.Normal(:,2), ptCloud2D.Normal(:,1)];
ptCloud2D.Location = ptCloud2D.Location + e2 .* t;
ptCloud2D.Normal = ptCloud2D.Normal + e3;
ptCloud2D.Normal = ptCloud2D.Normal./(ptCloud2D.Normal(:,1).^2 + ptCloud2D.Normal(:,2).^2);
pcNoisy = ptCloud2D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudDown = pcdownsample2D(ptCloud2D, percentage)
%pcdownsample2D Downsample a 2-D point cloud
% returns a downsampled point cloud with random sampling and without
% replacement.
ind = randi(ptCloud2D.Count, 1, round(ptCloud2D.Count*percentage));
location = ptCloud2D.Location(ind,:);
normal = ptCloud2D.Normal(ind,:);
ptCloudDown = pointCloud2D(location, normal);
end