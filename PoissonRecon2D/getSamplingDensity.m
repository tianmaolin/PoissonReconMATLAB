function weights = getSamplingDensity(samples, tree)
% Maolin Tian, 2018
global valueTable
weights = zeros(samples.Count, 1);
for s1 = 1:samples.Count
  n1 = samples.tree_ind(s1);
  for s2 = cell2mat(tree.sample_ind(tree.ngbr{n1}))'
    n2 = samples.tree_ind(s2);
    d = tree.depth(n2);
    dp = samples.Location(s2, :) - samples.Location(s1, :);
    dx = 1 + round((dp(1)-valueTable{d}(1, 1))*2^(d + 2));
    dy = 1 + round((dp(2)-valueTable{d}(1, 1))*2^(d + 2));
    Len = size(valueTable{d}, 1);
    if dx <= 0 || dx > Len || dy <= 0 || dy > Len
      continue;
    end
    F_i = valueTable{d}(dx, 2) * valueTable{d}(dy, 2) / (tree.width(n2))^2;
    weights(s1) = weights(s1) + F_i;
  end
end
