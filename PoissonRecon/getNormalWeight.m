function normalWeights = getNormalWeight(samples, tree, locationWeights)
%getWeight weight_s = \sum_{s\in S} 1/s.w^3 b_{s.w}(x-s.x) b_{s.w}(y-s.y)
%b_{s.w}(z-s.z)

global valueTable
normalWeights = zeros(samples.Count, 3);
for s1 = 1:samples.Count
    n1 = samples.tree_ind(s1);
    weigthSum = 0;
    for s2 = cell2mat(tree.sample_ind(tree.ngbr{n1}))'
        n2 = samples.tree_ind(s2);
        d = tree.depth(n2);
        dp = samples.Location(s2,:) - samples.Location(s1,:);
        dx = 1+round((dp(1) - valueTable{d}(1,1)) * 2^(d+2));
        dy = 1+round((dp(2) - valueTable{d}(1,1)) * 2^(d+2));
        dz = 1+round((dp(3) - valueTable{d}(1,1)) * 2^(d+2));
        Len = size(valueTable{d},1);
        if dx <= 0 || dx > Len || dy <= 0 || dy > Len || dz <= 0 || dz > Len
            continue;
        end
        
        F_i = valueTable{d}(dx,2) * valueTable{d}(dy,2) * valueTable{d}(dz,2);
        F_i = F_i / (tree.width(n2))^3;
        weigthSum = weigthSum + F_i / locationWeights(s2);
        F_i = F_i / locationWeights(s2) * samples.Normal(s2, :);
        normalWeights(s1, :) = normalWeights(s1, :) + F_i;
    end
    normalWeights(s1, 1) = norm(normalWeights(s1, :)) / weigthSum;
end
normalWeights = normalWeights(:, 1);