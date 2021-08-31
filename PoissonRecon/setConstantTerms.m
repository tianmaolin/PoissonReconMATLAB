function b = setConstantTerms(tree, samples, weight)
%setConstantTerms b_i = int V * d F_i dp
%
% Maolin Tian, 2018
global dotdTable

N = tree.Count;
int_F_dxF = cell(tree.maxDepth, tree.minDepth);
int_F_dyF = cell(tree.maxDepth, tree.minDepth);
int_F_dzF = cell(tree.maxDepth, tree.minDepth);
for d1 = tree.minDepth:tree.maxDepth
  for d2 = tree.minDepth:tree.maxDepth
    dx = dotdTable{d1, d2}(1, 1):2^(-max(d1, d2) - 1):dotdTable{d1, d2}(end, 1);
    [dx, dy, dz] = meshgrid(dx, dx, dx);
    int_F_dxF{d1, d2} = reshape(int_Fs_dxFi(d1, d2, dx, dy, dz), size(dx));
    int_F_dyF{d1, d2} = reshape(int_Fs_dyFi(d1, d2, dx, dy, dz), size(dx));
    int_F_dzF{d1, d2} = reshape(int_Fs_dzFi(d1, d2, dx, dy, dz), size(dx));
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
    p = samples.Location(s, :);
    normal = samples.Normal(s, :);
    o = tree.center(n, :);
    if max(abs(p - o)) >= 1.5 * (w1 + w2)
      continue;
    end

    dx = p(1) - o(1);
    dx = 1 + round((dx-dotdTable{d1, d2}(1, 1))*2^(max(d1, d2) + 1));
    dy = p(2) - o(2);
    dy = 1 + round((dy-dotdTable{d1, d2}(1, 1))*2^(max(d1, d2) + 1));
    dz = p(3) - o(3);
    dz = 1 + round((dz-dotdTable{d1, d2}(1, 1))*2^(max(d1, d2) + 1));
    Len = size(dotdTable{d1, d2}, 1);
    if dx <= 0 || dx > Len || dy <= 0 || dy > Len || dz <= 0 || dz > Len
      continue;
    end
    temp = [int_F_dxF{d1, d2}(dy, dx, dz); int_F_dyF{d1, d2}(dy, dx, dz); int_F_dzF{d1, d2}(dy, dx, dz)];
    t_n = normal * temp;
    b(n) = b(n) + t_n / weight(s);
  end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_B_dxB = int_Fs_dxFi(d1, d2, dx, dy, dz)
%int_Fs_dxFi int F1 dxF2 dx
if d1 >= d2
  d1S = d1;
  d2S = d2;
else
  d1S = d2;
  d2S = d1;
end
global dotTable dotdTable

if abs(dx) >= -dotdTable{d1, d2}(1, 1) | abs(dy) >= -dotTable{d1S, d2S}(1, 1) ...
    | abs(dz) >= -dotTable{d1S, d2S}(1, 1)
  int_B_dxB = 0;
  return;
end
m1 = round(dx.*2^(d1S + 1)+size(dotdTable{d1, d2}, 1)/2);
m2 = round(dy.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
m3 = round(dz.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
int_B_dxB = dotdTable{d1, d2}(m1, 2) .* dotTable{d1S, d2S}(m2, 2) .* dotTable{d1S, d2S}(m3, 2);
int_B_dxB = int_B_dxB * 8^d1 * 8^d2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_B_dyB = int_Fs_dyFi(d1, d2, dx, dy, dz)
%int_Fs_dyFi int F1 dyF2 dx
if d1 >= d2
  d1S = d1;
  d2S = d2;
else
  d1S = d2;
  d2S = d1;
end
global dotTable dotdTable

if abs(dx) >= -dotTable{d1S, d2S}(1, 1) | abs(dy) >= -dotdTable{d1, d2}(1, 1) ...
    | abs(dz) >= -dotTable{d1S, d2S}(1, 1)
  int_B_dyB = 0;
  return;
end
m1 = round(dx.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
m2 = round(dy.*2^(d1S + 1)+size(dotdTable{d1, d2}, 1)/2);
m3 = round(dz.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
int_B_dyB = dotTable{d1S, d2S}(m1, 2) .* dotdTable{d1, d2}(m2, 2) .* dotTable{d1S, d2S}(m3, 2);
int_B_dyB = int_B_dyB * 8^d1 * 8^d2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function int_B_dzB = int_Fs_dzFi(d1, d2, dx, dy, dz)
%int_Fs_dzFi int F1 dzF2 dx
if d1 >= d2
  d1S = d1;
  d2S = d2;
else
  d1S = d2;
  d2S = d1;
end
global dotTable dotdTable

if abs(dx) >= -dotTable{d1S, d2S}(1, 1) | abs(dy) >= -dotTable{d1S, d2S}(1, 1) ...
    | abs(dz) >= -dotdTable{d1, d2}(1, 1)
  int_B_dzB = 0;
  return;
end
m1 = round(dx.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
m2 = round(dy.*2^(d1S + 1)+size(dotTable{d1S, d2S}, 1)/2);
m3 = round(dz.*2^(d1S + 1)+size(dotdTable{d1, d2}, 1)/2);
int_B_dzB = dotTable{d1S, d2S}(m1, 2) .* dotTable{d1S, d2S}(m2, 2) .* dotdTable{d1, d2}(m3, 2);
int_B_dzB = int_B_dzB * 8^d1 * 8^d2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
