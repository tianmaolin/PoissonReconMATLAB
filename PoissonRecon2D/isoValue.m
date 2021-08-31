function v = isoValue(tree, samples, x)
% Maolin Tian, 2018
global valueTable
v = 0;
for s = 1:samples.Count
  n = samples.tree_ind(s);
  for m = tree.ngbr{n}'
    d = tree.depth(m);
    dp = samples.Location(s, :) - tree.center(m, :);
    dx = 1 + round((dp(1)-valueTable{d}(1, 1))*2^(d + 2));
    dy = 1 + round((dp(2)-valueTable{d}(1, 1))*2^(d + 2));
    Len = size(valueTable{d}, 1);
    if dx <= 0 || dx > Len || dy <= 0 || dy > Len
      continue;
    end
    F_i = valueTable{d}(dx, 2) * valueTable{d}(dy, 2) / (tree.width(m))^2;
    v = v + x(m) * F_i;
  end
end
v = v / samples.Count;
