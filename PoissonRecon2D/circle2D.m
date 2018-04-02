function circle = circle2D(N)
%circle2D make 2D pointCloud of a circle

% % circle
% dt = pi / round(N);
% alpha = [0:3*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:dt:2*pi-dt];
% % alpha = 0:2*dt:2*pi-dt;
% x = cos(alpha);
% y = sin(alpha);
% location = [x(:), y(:)];
% normal = - [x(:), y(:)];
% circle = pointCloud2D(location, normal);

% triangle
% x0 = 0:1/N:1;
% y0 = 1 - x0;
% x = [x0, x0, zeros(1, length(x0))];
% y = [zeros(1, length(x0)), y0, y0];
% location = [x(:), y(:)];
% n1 = repmat([0,1], length(x0),1);
% n2 = repmat([-sqrt(2)/2, -sqrt(2)/2],length(x0),1);
% n3 = repmat([1,0], length(x0),1);
% normal = [n1; n2; n3];
% circle = pointCloud2D(location, normal);

% circle with error
dt = pi / round(N);
alpha = [0:3*dt:pi/2-dt, pi/2:2*dt:3/2*pi-dt, 3/2*pi:10*dt:2*pi-dt];
% alpha = 0:2*dt:2*pi-dt;
x = cos(alpha);
y = sin(alpha);
location = [x(:), y(:)];
loc_error = 0.1*rand(size(location));
normal = - [x(:), y(:)];
nor_error = 0.1*rand(size(normal));
circle = pointCloud2D(location+loc_error, normal+nor_error);
