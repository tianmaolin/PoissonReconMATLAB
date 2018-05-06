function weights = getWeight(samples, minDepth, maxDepth)
%getWeight weight_s = \sum_{s\in S} 1/s.w^3 b_{s.w}(x-s.x) b_{s.w}(y-s.y)
%b_{s.w}(z-s.z)

if maxDepth < minDepth
    disp('maxDepth < minDepth !')
    return;
end

global valueTable
[tree,samples] = setTree(samples, minDepth, maxDepth);
weights = zeros(samples.Count, 1);
for s1 = 1:samples.Count
    n1 = samples.tree_ind(s1);
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
        weights(s1) = weights(s1) + F_i;
    end
end
