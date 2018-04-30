function A = getCoefficients(grid, B_dim)
%getCoefficients Coefficients_{ij} = int d B_i * d B_j dp
% Dirichlet boundary condition = 0
%
% Maolin Tian, Tongji University, 2018

N = 1/grid.width;
if N ~= floor(N)
    error('grid.width must = 1/N');
end

A_ij = zeros(5,5);
for i = 1:5
    for j = 1:5
        A_ij(i,j) = int_dB_dB([i-3,j-3], B_dim) / grid.width^4;
    end
end

% map = @(i,j) i + (j - 1) * N;
% map_i = @(n) mod(n,N);
% map_j = @(n) ceil(n/N);
I = zeros(125 * N ^ 2, 1);
J = zeros(125 * N ^ 2, 1);
V = zeros(125 * N ^ 2, 1);

m = 1:N^2;
% 0 Dirichlet, A_ii = 1, b = 0
boundm = m(mod(m,N) == 1 | mod(m,N) == N | ceil(m/N) == 1 | ceil(m/N) == N);
for j = 1:5
    for i = 1:5
        indT = (i-1 + (j-1) * 5) * N^2;
        I(m + indT) = m;
        J(m + indT) = m + (i-3) + (j-3) * N;
        V(m + indT) = repmat(A_ij(i,j), 1, N^2);
        I(boundm + indT) = [];
        J(boundm + indT) = [];
        V(boundm + indT) = [];
    end
end
I(boundm) = boundm;
J(boundm) = boundm;
V(boundm) = 1;
nanJ = J <= 0 | J > N * N;
I(nanJ) = [];
J(nanJ) = [];
V(nanJ) = [];

% m = 1;
% for i = 1:N
%     for j = 1:N
%         if i == 1 || i == N || j == 1 || j == N
%             % 0 Dirichlet, A_ii = 1, b = 0
%             I(m) = map(i,j);
%             J(m) = map(i,j);
%             V(m) = 1;
%             m = m + 1;
%             continue;
%         end
%         for a = 1:5
%             for b = 1:5
%                 if i+a-3 <= 0 || i+a-3 > N || j+b-3 <= 0 || j+b-3 > N
%                     continue;
%                 end
%                 if A_ij(a,b) == 0
%                     continue;
%                 end
%                 I(m) = map(i,j);
%                 J(m) = map(i+a-3,j+b-3);
%                 V(m) = A_ij(a,b);
%                 m = m + 1;
%             end
%         end
%     end
% end
% i = find(I == 0, 1);
% I = I(1 : i - 1);
% J = J(1 : i - 1);
% V = V(1 : i - 1);

A = sparse(I, J, V);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A_ij = int_dB_dB(trans, B_dim)
%int_dB_dB int dB * dB dx

% basic function
if B_dim == 2
    B = bspline([-1.5, -0.5, 0.5, 1.5]);
elseif B_dim == 1
	B = bspline([-1, 0, 1]);
else
    error('B_dim should be equal to 1 or 2.')
end
dB = fnder(B);
int_dxB_dxB = fn_int_F_Ft(dB, dB, trans(1)) * fn_int_F_Ft( B,  B, trans(2));
int_dyB_dyB = fn_int_F_Ft( B,  B, trans(1)) * fn_int_F_Ft(dB, dB, trans(2));
A_ij = int_dxB_dxB + int_dyB_dyB;
end
