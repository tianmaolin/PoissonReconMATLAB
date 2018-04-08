function v = getIsoValue(grid, samples, weight, x, B_dim, B_div)
%getLinearBsplineSum v = {sum_i X(s_i)/w(s_i)} / {sum_i 1/w(s_i)}
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

X = zeros(samples.Count, 1);
for s = 1 : samples.Count
    p = samples.Location(s, :);
    X(s) = getLinearBsplineSumP(grid, x, p, B_dim, B_div, y);
end

v = sum(X ./ weight) / sum(1./weight);

end



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