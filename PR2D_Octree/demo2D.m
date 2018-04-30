% Shows a couple of sample reconstruction using PR - Poisson Surface
% Reconstruction
%
% Maolin Tian, Tongji University, 2018

% ptCloud = ptCloudExample2D('circle', 1000);
% ptCloud = ptCloudExample2D('Armadillo');
ptCloud = ptCloudExample2D('Dragon');
% ptCloud = ptCloudDownsaple;
depth = 5; % grid.size = [2^depth, 2^depth]
verbose = false;

figure
plot(ptCloud.Location(:,1), ptCloud.Location(:,2), '.')
title('Input Point Cloud')

Contour = poissonRecon2D(ptCloud, depth, verbose);
if Contour(2,1) == size(Contour, 2) - 1
    disp('Number of Isoline pieces: ')
else
    disp('Number of Main Isoline pieces: ')
end
disp(Contour(2, 1))

% figure, hold on
% plot(Contour(1,2:end),Contour(2,2:end))
% dt = pi / 100;
% alpha = 0:2*dt:2*pi;
% xcir = cos(alpha);
% ycir = sin(alpha);
% plot(xcir, ycir)
% legend('recon','circle')
