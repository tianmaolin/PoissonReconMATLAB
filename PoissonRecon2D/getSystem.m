function A = getSystem(w)
% Dirichlet boundary condition = 0

N = 1/w;
if N ~= floor(N)
    disp('1/w is not int')
    return
end
A_ii = zeros(5,5);
for i = 1:5
    for j = 1:5
        A_ii(i,j) = bformInt([i-3,j-3],w);
    end
end

map = @(i,j) i + (j-1)*N;
% inv_map_i = @(n) mod(n,N);
% inv_map_j = @(n) (mod(n,N*N) - mod(n,N))/N + 1;
% inv_map_k = @(n) (n - mod(n,N*N))/N/N + 1;
I = zeros(125*N^2,1);
J = zeros(125*N^2,1);
V = zeros(125*N^2,1);
u = 1;
% TODO: use kron() instead of 4 for...end
for i = 1:N
    for j = 1:N 
        if i == 1 || i == N || j == 1 || j == N
            % 0 Dirichlet, A_ii = 1, b = 0
            I(u) = map(i,j);
            J(u) = map(i,j);
            V(u) = 1;
            u = u + 1;
        else
            for a = 1:5
                for b = 1:5
                    if i+a-3 <= 0 || i+a-3 > N || j+b-3 <= 0 || j+b-3 > N
                        continue;
                    end
                    if A_ii(a,b) == 0
                        continue;
                    end
                    I(u) = map(i,j);
                    J(u) = map(i+a-3,j+b-3);
                    V(u) = A_ii(a,b);
                    u = u + 1;
                end
            end
        end
    end
end
i = find(I==0, 1);
I = I(1:i-1);
J = J(1:i-1);
V = V(1:i-1);
A = sparse(I, J, V);
end

function A_ii = bformInt(trans,w)
A_ii = 0;
% F = bspline([-1.5,-0.5,0.5,1.5]); % basic function
F = bspline([-1,0,1]); % basic function
dF = fnder(F);
trans = - trans;

A_dx = 1;
t = trans(1);
F_int = fnmult(dF,fntrans(dF,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_dx = A_dx * F_int;
t = trans(2);
F_int = fnmult(F,fntrans(F,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_dx = A_dx * F_int;
A_dx = A_dx/w^4;
A_ii = A_ii + A_dx;

A_dy = 1;
t = trans(1);
F_int = fnmult(F,fntrans(F,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_dy = A_dy * F_int;
t = trans(2);
F_int = fnmult(dF,fntrans(dF,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_dy = A_dy * F_int;
A_dy = A_dy/w^4;
A_ii = A_ii + A_dy;
end