% show quadtree
% tree is in poissonRecon2D()
figure, hold on
for i = 1:tree.Count
    rectangle('Position',[tree.center(i,1)-tree.width(i)/2 ...
        tree.center(i,2)-tree.width(i)/2 tree.width(i) tree.width(i)])
end
plot(ptCloud2d.Location(:,1),ptCloud2d.Location(:,2),'b.')
axis equal

% % show b and x of 2-D grid data
% N = 1/w;
% figure
% U1 = w/2:w:1+2*w;
% [U1,V1]= meshgrid(U1, U1);
% U1 = double((U1 - 0.5) * scale - T(1));
% V1 = double((V1 - 0.5) * scale - T(2));
% Z1 = reshape(X,N,N);
% Z1 = [Z1; zeros(2,size(Z1,2))];
% Z1 = [Z1, zeros(size(Z1,1),2)];
% Z1(end,end) = -max(X);
% surf(U1, V1, Z1)
% shading interp
% title('\chi')
% axis([-1.477, 1.477, -1.477, 1.477])
% figure
% Z1 = reshape(b,N,N);
% surf(U, V, Z1)
% shading interp
% title('b')
% axis([-1.477, 1.477, -1.477, 1.477])
