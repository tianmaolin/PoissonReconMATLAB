function v = isoValue(tree, samples, x)

global valueTable
v = 0;
for s = 1:samples.Count
    n = samples.tree_ind(s);
    for m = tree.ngbr{n}'
        d = tree.depth(m);
        dp = samples.Location(s,:) - tree.center(m,:);
        dx = 1+round((dp(1) - valueTable{d}(1,1)) * 2^(d+2));
        dy = 1+round((dp(2) - valueTable{d}(1,1)) * 2^(d+2));
        dz = 1+round((dp(3) - valueTable{d}(1,1)) * 2^(d+2));
        Len = size(valueTable{d},1);
        if dx <= 0 || dx > Len || dy <= 0 || dy > Len || dz <= 0 || dz > Len
            continue;
        end
        F_i = valueTable{d}(dx,2) * valueTable{d}(dy,2) * valueTable{d}(dz,2);
        F_i = F_i * 8^(tree.depth(m));
        v = v + x(m) * F_i;
    end
end
v = v / samples.Count;
