function b = setConstantTerms(tree, samples, weight)
%setConstantTerms b_i = int V * d F_i dp
%
% Maolin Tian, Tongji University, 2018
global dotdTable

N = tree.Count;
int_F_dxF = cell(tree.maxDepth,tree.minDepth);
int_F_dyF = cell(tree.maxDepth,tree.minDepth);
for d1 = tree.minDepth:tree.maxDepth
    for d2 = tree.minDepth:tree.maxDepth
        dx = dotdTable{d1,d2}(1,1):2^(-max(d1,d2)-1):dotdTable{d1,d2}(end,1);
        [dx,dy] = meshgrid(dx,dx);
        int_F_dxF{d1,d2} = reshape(int_Fs_dxFi(d1, d2, dx, dy),size(dx));
        int_F_dyF{d1,d2} = reshape(int_Fs_dyFi(d1, d2, dx, dy),size(dx));
    end
end

b = zeros(N, 1);
for n = 1:N
    w1 = tree.width(n);
    for s = cell2mat(tree.sample_ind(tree.ngbr{n}))'
        m = samples.tree_ind(s);
        d1 = tree.depth(n);
        d2 = tree.depth(m);
        w2 = tree.width(m);
        p = samples.Location(s,:);
        normal = samples.Normal(s,:);
        o = tree.center(n,:);
        if max(abs(p - o)) >= 1.5 * (w1 + w2)
            continue;
        end
        
        dx = p(1) - o(1);
        dx = 1+round((dx - dotdTable{d1,d2}(1,1)) * 2^(max(d1,d2)+1));
        dy = p(2) - o(2);
        dy = 1+round((dy - dotdTable{d1,d2}(1,1)) * 2^(max(d1,d2)+1));
        Len = size(dotdTable{d1,d2},1);
        if dx <= 0 || dx > Len || dy <= 0 || dy > Len
            continue;
        end
        temp = [int_F_dxF{d1,d2}(dy, dx); int_F_dyF{d1,d2}(dy, dx)];
        t_n = normal * temp;
        b(n) = b(n) + t_n / weight(s);
    end
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_B_dxB = int_Fs_dxFi(d1, d2, dx, dy)
%int_dF_dF int F1 dxF2 dx
if d1 >= d2
    d1S = d1;
    d2S = d2;
else
    d1S = d2;
    d2S = d1;
end
global dotTable dotdTable

if abs(dx) >= - dotdTable{d1,d2}(1,1) | abs(dy) >= - dotTable{d1S,d2S}(1,1)
    int_B_dxB = 0;
else
    m = round(dx .* 2^(d1S+1) + size(dotdTable{d1,d2},1)/2);
    n = round(dy .* 2^(d1S+1) + size(dotTable{d1S,d2S},1)/2);
    int_B_dxB = dotdTable{d1,d2}(m,2) .* dotTable{d1S,d2S}(n,2);
end
int_B_dxB = int_B_dxB * 4^d1 * 4^d2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_B_dyB = int_Fs_dyFi(d1, d2, dx, dy)
%int_dF_dF int F1 dyF2 dx
if d1 >= d2
    d1S = d1;
    d2S = d2;
else
    d1S = d2;
    d2S = d1;
end
global dotTable dotdTable

if abs(dx) >= - dotTable{d1S,d2S}(1,1) | abs(dy) >= - dotdTable{d1,d2}(1,1)
    int_B_dyB = 0;
else
    m = round(dx .* 2^(d1S+1) + size(dotTable{d1S,d2S},1)/2);
    n = round(dy .* 2^(d1S+1) + size(dotdTable{d1,d2},1)/2);
    int_B_dyB = dotTable{d1S,d2S}(m,2) .* dotdTable{d1,d2}(n,2);
end
int_B_dyB = int_B_dyB * 4^d1 * 4^d2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
