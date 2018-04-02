% %
% w = 0.1;
% c = [0.05, 0.05, 0.05; 0.15, 0.05, 0.05; 0.25, 0.05, 0.05; 0.35, 0.05, 0.05];
% F = @(x)(x >= -1.5) .* (x < -0.5) .* 1/2 .* (x + 3/2).^2 + ...
%         (x >= -0.5) .* (x <  0.5) .* (-x.^2 + 3/4) + ...
%         (x >=  0.5) .* (x <  1.5) .* 1/2 .* (x - 3/2) .^2;
% % B = @(k,x,y,z) F((x-c(k,1))/w) .* F((y-c(k,2))/w) .* F((z-c(k,3))/w) / (w^3);
% dF = @(x) (x >= -1.5) .* (x < -0.5) .* (x + 3/2) +...
%           (x >= -0.5) .* (x <  0.5) .* (-2 * x) +...
%           (x >=  0.5) .* (x <  1.5) .* (x - 3/2);
% dxB = @(k,x,y,z) dF((x-c(k,1))/w) .*  F((y-c(k,2))/w) .*  F((z-c(k,3))/w) / (w^4);
% dyB = @(k,x,y,z)  F((x-c(k,1))/w) .* dF((y-c(k,2))/w) .*  F((z-c(k,3))/w) / (w^4);
% dzB = @(k,x,y,z)  F((x-c(k,1))/w) .*  F((y-c(k,2))/w) .* dF((z-c(k,3))/w) / (w^4);
% A = zeros(3,3);
% for k = 1:4
%     fun = @(x,y,z) dxB(1,x,y,z) .* dxB(k,x,y,z)...
%         + dyB(1,x,y,z) .* dyB(k,x,y,z)+ dzB(1,x,y,z) .* dzB(k,x,y,z);
%     a = c(1)-1.5*w;
%     b = c(1)+1.5*w;
%     A(k) = integral3(fun,a,b,a,b,a,b,'AbsTol',1e-5,'RelTol',1e-1);
% end
% % fun = @(x,y,z) dxB(1,x,y,z) .* dxB(k,x,y,z)+ dyB(1,x,y,z) .* dyB(k,x,y,z)+ dzB(1,x,y,z) .* dzB(k,x,y,z);
% % fun = @(x,y) fun(x,y,0.05);
% % [x,y] = meshgrid(-0.1:0.01:0.2,-0.1:0.01:0.2);
% % z = fun(x,y);
% % mesh(x,y,z)


function A = getSystem(w)
% Dirichlet boundary condition = 0

N = 1/w;
if N ~= floor(N)
    disp('1/w is not int')
    return
end
A_ii = zeros(5,5,5);
% A(3,3,3) = bformInt([0,0,0],w);
% A(4,3,3) = bformInt([1,0,0],w);
% A(5,3,3) = bformInt([2,0,0],w);
% A(4,4,3) = bformInt([1,1,0],w);
% A(5,4,3) = bformInt([2,1,0],w);
% A(5,5,3) = bformInt([2,2,0],w);
% 
% A(4,4,4) = bformInt([1,1,1],w);
% A(5,4,4) = bformInt([2,1,1],w);
% A(5,5,4) = bformInt([2,2,1],w);
% 
% A(4,4,5) = bformInt([1,1,2],w);
% A(5,4,5) = bformInt([2,1,2],w);
% A(5,5,5) = bformInt([2,2,2],w);
for i = 1:5
    for j = 1:5
        for k = 1:5
            A_ii(i,j,k) = bformInt([i-3,j-3,k-3],w);
        end
    end
end

map = @(i,j,k) i + (j-1)*N + (k-1)*N*N;
% inv_map_i = @(n) mod(n,N);
% inv_map_j = @(n) (mod(n,N*N) - mod(n,N))/N + 1;
% inv_map_k = @(n) (n - mod(n,N*N))/N/N + 1;
I = zeros(125*N^3,1);
J = zeros(125*N^3,1);
V = zeros(125*N^3,1);
u = 1;
for i = 1:N
    for j = 1:N
        for k = 1:N % TODO: use kron() instead of 6 for...end
            if i == 1 || i == N || j == 1 || j == N || k == 1 || k == N
                % 0 Dirichlet, A_ii = 1, b = 0
                I(u) = map(i,j,k);
                J(u) = map(i,j,k);
                V(u) = 1;
                u = u + 1;
            else
                for a = 1:5
                    for b = 1:5
                        for c = 1:5
                            if i+a-3 > 0 && i+a-3 <= N && j+b-3 > 0 && j+b-3 <= N && k+c-3 > 0 && k+c-3 <= N
                                if A_ii(a,b,c) == 0
                                    continue;
                                end
                                I(u) = map(i,j,k);
                                J(u) = map(i+a-3,j+b-3,k+c-3);
                                V(u) = A_ii(a,b,c);
                                u = u + 1;
                            end
                        end
                    end
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
F = bspline([-1.5,-0.5,0.5,1.5]); % basic function
% F = bspline([-1,0,1]); % basic function
dF = fnder(F);
trans = - trans;

A_i = 1;
t = trans(1);
F_int = fnmult(dF,fntrans(dF,t));
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
A_i = A_i/w^5;
A_ii = A_ii + A_i;

A_i = 1;
t = trans(1);
F_int = fnmult(F,fntrans(F,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_i = A_i * F_int;
t = trans(2);
F_int = fnmult(dF,fntrans(dF,t));
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
A_i = A_i/w^5;
A_ii = A_ii + A_i;

A_i = 1;
t = trans(1);
F_int = fnmult(F,fntrans(F,t));
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
F_int = fnmult(dF,fntrans(dF,t));
if F_int.pieces == 0
    F_int = 0;
else
    F_int = fnval(fnint(F_int),F_int.breaks(end));
end
A_i = A_i * F_int;
A_i = A_i/w^5;
A_ii = A_ii + A_i;
end