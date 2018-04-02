function b = getb(grid, weight, div)
%getb get b

N = grid.size(1);
% map = @(i,j,k) i + (j-1)*N;
map_ij = @(n) [rem(n, N),floor(n / N) + 1];
w = grid.width;

dx = 1 / div;
x = [-3:dx:-dx, 0:dx:3];
ppformInt_c = ceil(length(x) / 2);
[X, Y] = meshgrid(x,x);
ppformInt = zeros(length(x), length(x));
for a = 1:length(x)
    for b = 1:length(x)
        ppformInt(a, b) = bformInt([X(a,b), Y(a,b)], w, div);
    end
end

grid_neighbor = cell(N^2,1);
for n = 1:N^2
    o = (map_ij(n) - 0.5) * grid.width;
    a = find(abs(grid.Location(:,1) - o(1))<= 3 * w);
    b = find(abs(grid.Location(:,2) - o(2))<= 3 * w);
    grid_neighbor{n} = intersect(a,b)';
end

b = zeros(N^2, 1);
for n = 1:N^2
    if isempty(grid_neighbor{n})
        continue;
    end
    for s = grid_neighbor{n}
        p = grid.Location(s,:);
        normal = grid.Normal(s,:);
        o = (map_ij(n) - 0.5) * grid.width;
        if max(abs(p - o)) >= 3 * w
            continue;
        end
        t_x = (p(1)-o(1))/w; 
        t_y = (p(2)-o(2))/w;
        t_m = round(t_y/dx) + ppformInt_c;
        t_n = round(t_x/dx) + ppformInt_c;
        temp = [ppformInt(t_m,t_n); ppformInt(t_n,t_m)];
        t_n = normal * temp;
        b(n) = b(n) + t_n / weight(s);
    end
end
end

function A_i = bformInt(trans,w,div)
% F = bspline([-1.5,-0.5,0.5,1.5]); % basic function
F = bspline([-1,0,1]); % basic function
F = fndiv(F, div);
dF = fnder(F);
trans = - trans;

A_i = 1;
t = trans(1);
F_int = fnmult(F,fntrans(dF,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_i = A_i * F_int;

t = trans(2);
F_int = fnmult(F,fntrans(F,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_i = A_i * F_int;
A_i = A_i/w^3;
end
