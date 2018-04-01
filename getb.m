function b = getb(grid, weight, div)
%getb get b

N = grid.size(1);
% map = @(i,j,k) i + (j-1)*N + (k-1)*N*N;
map_ijk = @(n) [rem(n, N),floor(rem(n, N*N) / N),floor(n / N / N)];
w = grid.width;

dx = 1 / div;
x = [-3:dx:-dx, 0:dx:3];
ppformInt_c = ceil(length(x) / 2);
[X, Y, Z] = meshgrid(x,x,x);
ppformInt = zeros(length(x), length(x), length(x));
for a = 1:length(x)
    for b = 1:length(x)
        for c = 1:length(x)
            ppformInt(a, b, c) = bformInt([X(a,b,c), Y(a,b,c), Z(a,b,c)], w, div);
        end
    end
end

grid_neighbor = cell(N^3,1);
for n = 1:N^3
    o = (map_ijk(n) - 0.5) * grid.width;
    a = find(abs(grid.Location(:,1) - o(1))<= 3 * w);
    b = find(abs(grid.Location(:,2) - o(2))<= 3 * w);
    c = find(abs(grid.Location(:,3) - o(3))<= 3 * w);
    a = intersect(a,b);
    grid_neighbor{n} = intersect(a,c)';
end

b = zeros(N^3, 1);
for n = 1:N^3
    if isempty(grid_neighbor{n})
        continue;
    end
    for s = grid_neighbor{n}
        p = grid.Location(s,:);
        normal = grid.Normal(s,:);
        o = (map_ijk(n) - 0.5) * grid.width;
        if max(abs(p - o)) >= 3 * w
            continue;
        end
        temp = [ppformInt(round((p(1)-o(1))/w/dx) + ppformInt_c, round((p(2)-o(2))/w/dx) + ppformInt_c, round((p(3)-o(3))/w/dx) + ppformInt_c);
                ppformInt(round((p(2)-o(2))/w/dx) + ppformInt_c, round((p(1)-o(1))/w/dx) + ppformInt_c, round((p(3)-o(3))/w/dx) + ppformInt_c);
                ppformInt(round((p(3)-o(3))/w/dx) + ppformInt_c, round((p(2)-o(2))/w/dx) + ppformInt_c, round((p(1)-o(1))/w/dx) + ppformInt_c)];
        b(n) = b(n) + normal * temp / weight(s,:);
    end
end
end

function A_i = bformInt(trans,w,div)
F = bspline([-1.5,-0.5,0.5,1.5]); % basic function
% F = bspline([-1,0,1]); % basic function
F = fndiv(F, div);
dF = fnder(F);

A_i = 1;
t = trans(1);
% TODO:fnmult()
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
t = trans(3);
F_int = fnmult(F,fntrans(F,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_i = A_i * F_int;
A_i = A_i/w^4;
end