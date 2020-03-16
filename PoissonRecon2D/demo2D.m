% Show some curve reconstructions using Poisson Surface Reconstruction
%
% Maolin Tian, 2018

% Input
% ptCloud = ptCloudExample2D('triangle', 2000,0.002);
% ptCloud = ptCloudExample2D('circle', 1000);
% ptCloud = ptCloudExample2D('Armadillo');
ptCloud = ptCloudExample2D('Dragon');
minDepth = 4; % grid size = [2^depth, 2^depth]
maxDepth = 8; % It is better to less than 10.
verbose = false;

ptCloud = pcNormalized(ptCloud);
ptCloudTrain = pointCloud2D(ptCloud.Location(1:2:end,:), ptCloud.Normal(1:2:end,:));
ptCloudTest = pointCloud2D(ptCloud.Location(2:2:end,:), ptCloud.Normal(2:2:end,:));
figure
plot(ptCloudTrain.Location(:,1), ptCloudTrain.Location(:,2), '.')
title('Input Point Cloud')
% disp(['Depth:  ',  num2str(maxDepth)])

% Poisson Surface Reconstruction 2D
Contour = poissonRecon2D(ptCloudTrain, minDepth, maxDepth, verbose);
% if Contour(2,1) == size(Contour, 2) - 1
%     disp('Number of Isoline Pieces: ')
% else
%     disp('Number of Main Isoline Pieces: ')
% end
% disp(Contour(2, 1))

% % Error
% L = 1;
% V = [];
% figure, hold on
% while(true)
%     plot(Contour(1,L + 1:L + Contour(2,L)),Contour(2,L + 1:L + Contour(2,L)),'g')
%     V = [V, Contour(:,L + 1:L + Contour(2,L))];
%     L = L + Contour(2,L) + 1;
%     if L >= size(Contour,2)
%         break
%     end
% end
% [error, dist] = getError(V', ptCloudTest.Location);
% hold on, axis equal
% % plot(V(:,1), V(:,2), '.')
% plot(ptCloudTest.Location(dist>1.2*error,1), ptCloudTest.Location(dist>1.2*error,2), '*')
% title('The Points Contributing to Main Error')
% figure
% hist(dist)
% title('Error Distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ptCloud2D = ptCloudExample2D(model, N, noise)
%ptCloud2D generate 2D pointCloud examples
%
% Maolin Tian, 2018

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
   
elseif strcmp(model,'circle')
    dt = pi / round(N);
%    alpha = [0:3*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:10*dt:2*pi-dt];
%    alpha = [0:dt:pi/2-dt, pi/2:4*dt:3/2*pi-dt, 3/2*pi:10*dt:13/8*pi];
    alpha = 0:2*dt:2*pi-dt;
    x = cos(alpha);
    y = sin(alpha);
    location = [x(:), y(:)];
    normal = - [x(:), y(:)];
    ptCloud2D = pointCloud2D(location, normal);
        
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
    XZ_i = find(abs(ptCloud3D.Location(:,3) + 0.0045) < 0.0005 ); % ...
%         & ptCloud3D.Location(:,1) < -0.015 ...
%         & ptCloud3D.Location(:,1) > -0.105 ...
%         & ptCloud3D.Location(:,2) < 0.19 ...
%         & ptCloud3D.Location(:,2) > 0.12);
    location = ptCloud3D.Location(XZ_i,:);
    location = [location(:, 1), location(:, 2)];
    normal = ptCloud3D.Normal(XZ_i,:);
    normal = [normal(:, 1), normal(:, 2)];
    normal = normal./(normal(:,1).^2 + normal(:,2).^2);
    ptCloud2D = pointCloud2D(location, normal);
%     ptCloud2D = pcdownsample2D(ptCloud2D, 0.05);
   
end

if nargin > 2
    ptCloud2D = addNoise(ptCloud2D, noise);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pcNoisy = addNoise(ptCloud2D, e)
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

function ptCloud2D = pcNormalized(ptCloud2D)
length = sqrt(ptCloud2D.Normal(:,1).^2 + ptCloud2D.Normal(:,2).^2);
if find(length==0, 1)
    error('normal == 0 !')
end
ptCloud2D.Normal = ptCloud2D.Normal ./ length;
end
