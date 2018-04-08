function weight = getWeight(grid, samples, depth, B_dim, div)
%getWeight Estimate the number of points in a sample's neighborhood
%
% Maolin Tian, Tongji University, 2018

% downsample grid
N = 2^depth;
if grid.depth ~= depth
    grid = struct('depth', depth, 'width', 1/N, 'size', [N, N]);    
    samples.grid_ind = ceil(samples.Location / grid.width);
    samples.grid_ind = samples.grid_ind(:, 1) + (samples.grid_ind(:, 2) - 1) * N;
    grid.sample_ind = cell(N * N, 1);
    for n = 1 : N * N
        grid.sample_ind{n} = find(samples.grid_ind == n);
    end
end

ngbrT = repmat(-2:1:2, 5, 1) + repmat((-2*N:N:2*N)', 1, 5);
samples.neighbor = cell(samples.Count,1);
for s = 1:samples.Count
    ngbr = ngbrT + samples.grid_ind(s);
    for m = (ngbr(:))'
        if m < 1 || m > N * N
            continue
        end
        if isempty(grid.sample_ind{m})
            continue
        end
        samples.neighbor{s} = [samples.neighbor{s}, (grid.sample_ind{m})'];
    end
end

% basic function
if B_dim == 2
    B = bspline([-1.5, -0.5, 0.5, 1.5]);
elseif B_dim == 1
	B = bspline([-1, 0, 1]);
else
    error('B_dim should be equal to 1 or 2.')
end
B_end = B.breaks(end);
dx = 1 / div;
x = 0:dx:B_end;
y = fnval(B, x);

weight = zeros(samples.Count,1);
w = grid.width;
for n = 1:samples.Count
    for s = samples.neighbor{n}
        p = samples.Location(n,:);
        q = samples.Location(s,:);
        if max(abs(p - q)) >= B_end * w
            continue;
        end
        weight(n,:) = weight(n,:) + ...
            y(round(abs((p(1) - q(1)) / w)/dx) + 1) * ...
            y(round(abs((p(2) - q(2)) / w)/dx) + 1) / w^2;
    end
end
end
