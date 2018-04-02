function weight = getWeight(grid, depth, dx)
%getWeight get weights to estimate the number of points in a
% neighborhood of a sample
% log_2(depth) <= log_2(grid.depth)

% downsample grid to depth
if (grid.depth == depth)
    grid_low = grid;
else
    grid_low.depth = depth;
    grid_low.width = 1 / 2^grid_low.depth;
    grid_low.size = [1 / grid_low.width, 1 / grid_low.width];
    grid_low.Count = grid.Count;
    grid_low.Location = grid.Location;
    grid_low.Normal = grid.Normal;
    grid_low.gridijk = ceil(grid.Location / grid_low.width);
    grid_low.gridLocation = grid_low.gridijk * grid_low.width;
    grid_low.quadrant = (grid_low.gridLocation - grid_low.Location < grid_low.width/2);
    grid_low.neighbor = cell(grid_low.Count,1);
    for n = 1:grid_low.Count
        a = find(abs(grid_low.gridijk(:,1) - grid_low.gridijk(n,1))<= 2);
        b = find(abs(grid_low.gridijk(:,2) - grid_low.gridijk(n,2))<= 2);
        grid_low.neighbor{n} = intersect(a,b);
    end
end


B = bspline([-1.5, -0.5, 0.5, 1.5]);
x = 0:dx:1.5;
y = fnval(B, x);


weight = zeros(grid_low.Count,1);
for n = 1:grid_low.Count
    for s = grid_low.neighbor{n}'
        p = grid_low.Location(n,:);
        q = grid_low.Location(s,:);
        w = grid_low.width;
        if max(abs(p - q)) > 1.5 * w
            continue;
        end
        weight(n,:) = weight(n,:) + ...
            y(round(abs((p(1) - q(1)) / w)/dx) + 1) * ...
            y(round(abs((p(2) - q(2)) / w)/dx) + 1) / w^2;
    end
end
end
