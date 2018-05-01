function A = setCoefficients(tree)
%getCoefficients Coefficients_{ij} = int d B_i * d B_j dp
% Dirichlet boundary condition = 0
%
% Maolin Tian, Tongji University, 2018

% TODO: when degree = 2,maxDepth = minDepth, why the result of Ax = 4 is odd?
global dotTable

N = tree.Count;

A_ij = cell(tree.maxDepth,tree.minDepth);
for d1 = tree.minDepth:tree.maxDepth
    for d2 = tree.minDepth:d1
        dx = dotTable{d1,d2}(1,1):2^(-d1-1):dotTable{d1,d2}(end,1);
        [dx,dy] = meshgrid(dx,dx);
        % TODO: need trans matrix for reshape?
        A_ij{d1,d2} = reshape(int_dF1_dF2(d1, d2, dx, dy),size(dx));
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
%         if tree.isbound(j)
%             continue;
%         end
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
            d1 = tree.depth(iSort);
            d2 = tree.depth(jSort);
        end
        Len = size(dotTable{d1,d2},1);
        dp = tree.center(iSort,:) - tree.center(jSort,:);
        dx1 = dp(1);
        dy1 = dp(2);
        dx1 = 1+round((dx1 - dotTable{d1,d2}(1,1)) * 2^(d1+1));
        dy1 = 1+round((dy1 - dotTable{d1,d2}(1,1)) * 2^(d1+1));
        if dx1 <= 0 || dx1 > Len || dy1 <= 0 || dy1 > Len
            continue;
        end
        V(con) = A_ij{d1,d2}(dy1, dx1);
    end
end
nanJ = J <= 0 | J > N * N;
I(nanJ) = [];
J(nanJ) = [];
V(nanJ) = [];
A = sparse(I, J, V);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A_ij = int_dF1_dF2(d1, d2, dx, dy)
%int_dF_dF int dF1 dF2 dx
if d1 < d2
    d_temp = d1;
    d1 = d2;
    d2 = d_temp;
    d_temp = dx;
    dx = dy;
    dy = d_temp;
end
global dotTable ddotdTable

if abs(dx) >= - ddotdTable{d1,d2}(1,1) | abs(dy) >= - dotTable{d1,d2}(1,1)
    int_dxB_dxB = 0;
else
    m = round(dx .* 2^(d1+1) + size(ddotdTable{d1,d2},1)/2);
    n = round(dy .* 2^(d1+1) + size(dotTable{d1,d2},1)/2);
    int_dxB_dxB = ddotdTable{d1,d2}(m,2) .* dotTable{d1,d2}(n,2);
end
if abs(dx) >= - dotTable{d1,d2}(1,1) | abs(dy) >= - ddotdTable{d1,d2}(1,1)
    int_dyB_dyB = 0;
else
    m = round(dx .* 2^(d1+1) + size(dotTable{d1,d2},1)/2);
    n = round(dy .* 2^(d1+1) + size(ddotdTable{d1,d2},1)/2);
    int_dyB_dyB = dotTable{d1,d2}(m,2) .* ddotdTable{d1,d2}(n,2);
end
A_ij = int_dxB_dxB + int_dyB_dyB;
A_ij = - A_ij * 4^d1 * 4^d2;
end
