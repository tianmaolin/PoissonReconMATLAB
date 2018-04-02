% N = 400;
% x = -1.5:3/N:1.5-1/N;
% [x,y] = meshgrid(x,x);
% z = zeros(size(x));
% for i = 1:N
%     for j = 1:N
%         r = x(i,j)^2+y(i,j)^2;
%         if r<=1.1 && r>=0.9
%             z(i,j) = 0.5*(sin(5*(r-0.9)*pi+pi/2)+1);
%         elseif r<0.9
%             z(i,j) = 1;
%         end
%     end
% end
% mesh(x,y,z)

% N = 64;
x = -1.5:3/N:1.5-1/N;
[x,y] = meshgrid(x,x);
z = zeros(size(x));
for i = 1:N
    for j = 1:N
        r = x(i,j)^2+y(i,j)^2;
        if r<=1.1 && r>=0.9
            z(i,j) = sin(5*(r-0.9)*pi+pi/2) * 50*pi*pi*r-cos(5*(r-0.9)*pi+pi/2)*10*pi;
        end
    end
end
b = z(:);