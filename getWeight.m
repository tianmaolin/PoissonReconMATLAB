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
    grid_low.size = [1 / grid_low.width, 1 / grid_low.width, 1 / grid_low.width];
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
        c = find(abs(grid_low.gridijk(:,3) - grid_low.gridijk(n,3))<= 2);
        a = intersect(a,b);
        grid_low.neighbor{n} = intersect(a,c);
    end
end


B = bspline([-1.5, -0.5, 0.5, 1.5]);
x = 0:dx:1.5;
y = fnval(B, x);


weight = zeros(grid_low.Count,1);
for n = 1:grid_low.Count
    for s = grid_low.neighbor{n,1}'
% %         if norm(grid_low.gridijk(s) - grid_low.gridijk(n), 'inf') > 3
% %             continue;
% %         end
%         % suppose the neighbor of s is all in gird (s is not on the boundary)
%         s_ijk = grid_low.gridijk(s,:);
%         o = [s_ijk; s_ijk + [1,0,0]; s_ijk + [0,1,0]; s_ijk + [1,1,0];...
%             s_ijk + [0,0,1]; s_ijk + [1,0,1]; s_ijk + [0,1,1]; s_ijk + [1,1,1],];
%         o = o + grid.quadrant(n,:) - 1;
%         os_c = grid_low.gridLocation(s,:) + grid_low.width * (grid_low.quadrant(s,:) - 1);
%         F8 = zeros(8,1);
%         for i = 1:8
%             w = 1 / 2^grid_low.depth;
%             p = grid_low.Location(n,:);
%             F8(i) = fnval(B, p(1)/w - o(i,1) + 0.5) * fnval(B, p(2)/w - o(i,2) + 0.5) * fnval(B, p(3)/w - o(i,3) + 0.5) / w^3;
%         end
%         w = zeros(8,1);
%         x = (grid_low.Location(s,1) - os_c(1))/grid_low.width + 0.5;
%         y = (grid_low.Location(s,2) - os_c(2))/grid_low.width + 0.5;
%         z = (grid_low.Location(s,3) - os_c(3))/grid_low.width + 0.5;
%         w(1) = (1 - x)*(1 - y)*(1 - z) ;
%         w(2) = x*(1 - y)*(1 - z) ;
%         w(3) = (1 - x)*y*(1 - z) ;
%         w(4) = x*y*(1 - z) ;
%         w(5) = (1 - x)*(1 - y)*z ;
%         w(6) = x*(1 - y)*z ;
%         w(7) = (1 - x)*y*z ;
%         w(8) = x*y*z ;
        p = grid_low.Location(n,:);
        q = grid_low.Location(s,:);
        w = grid_low.width;
        if max(abs(p - q)) > 1.5 * w
            continue;
        end
        weight(n,:) = weight(n,:) + ...
            y(round(abs((p(1) - q(1)) / w)/dx) + 1) * ...
            y(round(abs((p(2) - q(2)) / w)/dx) + 1) * ...
            y(round(abs((p(3) - q(3)) / w)/dx) + 1) / w^3;
    end
end
end
