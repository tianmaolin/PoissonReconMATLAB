function X = getLinearBsplineSum(grid, x, B_dim, grid_div, B_div)
%getLinearBsplineSum X = sum_i x_i * B_i
%
% Maolin Tian, Tongji University, 2018

% basic function
if B_dim == 2
    B = bspline([-1.5, -0.5, 0.5, 1.5]);
elseif B_dim == 1
	B = bspline([-1, 0, 1]);
else
    error('B_dim should be equal to 1 or 2.')
end
B_end = B.breaks(end);
dx = 1 / B_div;
y = fnval(B, 0:dx:B_end);

N = 2 ^ grid.depth;
M = grid_div * N;
% ngbrT = repmat(-2:1:2, 5, 1) + repmat((-2*N:N:2*N)', 1, 5);
% grid.neighbor = cell(N * N, 1);
% for n = 1 : N * N
%     ngbr = ngbrT + n;
%     for m = (ngbr(:))'
%         if m < 1 || m > N * N
%             continue
%         end
%         if isempty(grid.sample_ind{m})
%             continue
%         end
%         grid.neighbor{n} = [grid.neighbor{n}, (grid.sample_ind{m})'];
%     end
% end

map_ij = @(n) [rem(n, M),ceil(n / M)];
w = grid.width / grid_div;
X = zeros(M ^ 2, 1);
for m = 1 : M ^ 2
    p = map_ij(m) * w;
    X(m) = getLinearBsplineSumP(grid, x, p, B_dim, B_div, y);
%     n = ceil(m / grid_div / grid_div);
%     for s = grid.neighbor{n}
%         p = samples.Location(s,:);
%         q = map_ij(m) * w;
%         d = norm(p - q) / grid_div / w;
%         if d  > B_end
%             continue;
%         end
%         X(m) = X(m) + y(round(d / dx + 1)) * x(s);
%     end
end
end


% if floor(p_N / N) ~= p_N / N
%     Sum = 0;
%     return;
% end
% 
% w = 1/N;
% p = 0:1/p_N:1;
% [px, py] = meshgrid(p, p);
% px = px(:);
% py = py(:);
% map_xy = @(n) [rem(n, N),floor(n / N) + 1];
% 
% neighbor = cell(N^2,1);
% for n = 1:N^2
%     o = (map_xy(n) - 0.5) * w;
%     a = find(abs(px(:) - o(1)) < 2 * w);
%     b = find(abs(py(:) - o(2)) < 2 * w);
%     neighbor{n} = intersect(a,b)';
% end
% 
% % B = bspline([-1.5,-0.5,0.5,1.5]); % basic function
% B = bspline([-1,0,1]); % basic function
% F = zeros(B_N, B_N);
% for i = 1:B_N
%     for j = 1:B_N
%         Bx = (i-1) / B_N * 2;
%         By = (j-1) / B_N * 2;
%         if Bx > B.breaks(end) || By > B.breaks(end)
%             F(i,j) = 0;
%         else
%             F(i,j) = fnval(B, Bx) * fnval(B, By) / w^2;
%         end
%     end
% end
% 
% map_xy = @(n) [rem(n, p_N),floor(n / p_N) + 1];
% Sum = zeros(p_N ^ 2,1);
% for n = 1:p_N^2
%     for i = neighbor{n}
%         c = (map_xy(i) - 0.5) * w;
%         t = [abs(px(map_xy(n)) - c(1)) / w, abs(px(map_xy(n)) - c(2)) / w];
%         t = [floor(t(1) * p_N/2) + 1, floor(t(2) * p_N/2) + 1];
%         if max(t) > p_N
%             continue;
%         end
%         B = F(t(1)) * F(t(2)) / w^2;
%         Sum(n) = Sum(n) + x(n) * B;
%     end
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = getLinearBsplineSumP(grid, x, P, B_dim, B_div, y)
%getLinearBsplineSumP Get Linear B-spline Sum at point P
% X = sum_i x_i * B_i
N = 2 ^ grid.depth;
P_ind = ceil(P / grid.width);
P_ind = P_ind(1) + (P_ind(2) - 1) * N;
ngbrT = repmat(-2:1:2, 5, 1) + repmat((-2*N:N:2*N)', 1, 5);
ngbr = ngbrT + P_ind;
ngbr(ngbr < 1 | ngbr > N * N) = [];
ngbr = (ngbr(:))';
map_ij = @(n) [rem(n, N),ceil(n / N)];
w = grid.width;
X = 0;
if B_dim == 2
    B_end = 1.5;
else
    B_end = 1;
end
for n = ngbr
    o = (map_ij(n) - 0.5) * w;
    d = norm(P - o) / w;
    if d  > B_end
        continue;
    end
    X = X + y(round(d * B_div + 1)) * x(n);
end


end