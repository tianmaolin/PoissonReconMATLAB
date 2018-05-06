function A = setCoefficients(tree)
%setCoefficients Coefficients_{ij} = int d B_i * d B_j dp
% Dirichlet boundary condition = 0
%
% Maolin Tian, Tongji University, 2018

global dotTable

N = tree.Count;

A_ij = cell(tree.maxDepth,tree.minDepth);
for d1 = tree.minDepth:tree.maxDepth
    for d2 = tree.minDepth:d1
        dx = dotTable{d1,d2}(1,1):2^(-d1-1):dotTable{d1,d2}(end,1);
        [dx,dy,dz] = meshgrid(dx,dx,dx);
        A_ij{d1,d2} = reshape(int_dF1_dF2(d1, d2, dx, dy, dz),size(dx));
    end
end

I = zeros(100 * N, 1);
J = zeros(100 * N, 1);
V = zeros(100 * N, 1);
con = 0;
for i = 1:N
    if tree.isbound(i)
        con = con + 1;
        I(con) = i;
        J(con) = i;
        V(con) = 1;
        continue;
    end
    for j = tree.ngbr{i}'
        con = con + 1;
        I(con) = i;
        J(con) = j;
        iSort = i;
        jSort = j;
        d1 = tree.depth(i);
        d2 = tree.depth(j);
        if d1 < d2
            iSort = j;
            jSort = i;
            d_temp = d1;
            d1 = d2;
            d2 = d_temp;
        end
        Len = size(dotTable{d1,d2},1);
        dp = tree.center(iSort,:) - tree.center(jSort,:);
        dx1 = dp(1);
        dy1 = dp(2);
        dz1 = dp(3);
        dx1 = 1+round((dx1 - dotTable{d1,d2}(1,1)) * 2^(d1+1));
        dy1 = 1+round((dy1 - dotTable{d1,d2}(1,1)) * 2^(d1+1));
        dz1 = 1+round((dz1 - dotTable{d1,d2}(1,1)) * 2^(d1+1));
        if dx1 <= 0 || dx1 > Len || dy1 <= 0 || dy1 > Len || dz1 <= 0 || dz1 > Len
            continue;
        end
        V(con) = A_ij{d1,d2}(dy1, dx1, dz1);
    end
end
nanJ = J <= 0 | J > N * N;
I(nanJ) = [];
J(nanJ) = [];
V(nanJ) = [];
A = sparse(I, J, V);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A_ij = int_dF1_dF2(d1, d2, dx, dy, dz)
%int_dF_dF int dF1 dF2 dx
if d1 < d2
    d_temp = d1;
    d1 = d2;
    d2 = d_temp;
end
global dotTable ddotdTable

if abs(dx) >= - dotTable{d1,d2}(1,1) | abs(dy) >= - dotTable{d1,d2}(1,1)...
        | abs(dz) >= - dotTable{d1,d2}(1,1)
    A_ij = 0;
    return;
end

m1 = round(dx .* 2^(d1+1) + size(dotTable{d1,d2},1)/2);
m2 = round(dy .* 2^(d1+1) + size(dotTable{d1,d2},1)/2);
m3 = round(dz .* 2^(d1+1) + size(dotTable{d1,d2},1)/2);
int_dxB_dxB = ddotdTable{d1,d2}(m1,2) .*   dotTable{d1,d2}(m2,2) .*   dotTable{d1,d2}(m3,2);
int_dyB_dyB =   dotTable{d1,d2}(m1,2) .* ddotdTable{d1,d2}(m2,2) .*   dotTable{d1,d2}(m3,2);
int_dzB_dzB =   dotTable{d1,d2}(m1,2) .*   dotTable{d1,d2}(m2,2) .* ddotdTable{d1,d2}(m3,2);

A_ij = int_dxB_dxB + int_dyB_dyB + int_dzB_dzB;
A_ij = - A_ij * 8^d1 * 8^d2;
end
