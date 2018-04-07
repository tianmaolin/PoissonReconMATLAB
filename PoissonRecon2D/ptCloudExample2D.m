function ptCloud2D = ptCloudExample2D(N)
%ptCloud2D make 2D pointCloud example

% % circle
% dt = pi / round(N);
% alpha = [0:2*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:dt:2*pi-dt];
% % alpha = 0:2*dt:2*pi-dt;
% x = cos(alpha);
% y = sin(alpha);
% location = [x(:), y(:)];
% normal = - [x(:), y(:)];
% ptCloud2D = pointCloud2D(location, normal);

% % triangle
% x0 = 0:1/N:1;
% y0 = 1 - x0;
% x = [x0, x0, zeros(1, length(x0))];
% y = [zeros(1, length(x0)), y0, y0];
% location = [x(:), y(:)];
% n1 = repmat([0,1], length(x0),1);
% n2 = repmat([-sqrt(2)/2, -sqrt(2)/2],length(x0),1);
% n3 = repmat([1,0], length(x0),1);
% normal = [n1; n2; n3];
% ptCloud2D = pointCloud2D(location, normal);

% circle with error
dt = pi / round(N);
alpha = [0:3*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:10*dt:2*pi-dt];
% alpha = 0:2*dt:2*pi-dt;
x = cos(alpha);
y = sin(alpha);
location = [x(:), y(:)];
normal = - [x(:), y(:)];
ptCloud2D = pointCloud2D(location, normal);
ptCloud2D = AddNoise(ptCloud2D, 0.05);

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