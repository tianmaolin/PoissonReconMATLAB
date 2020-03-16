function curvatureField = getCurvatureField(samples, tree, locationWeights)
%getCurvatureField curvatureField(q) = sum_s F_s(q) / s.samplingDensity * s.N
% / sum_s F_s(q) / s.locationWeights
% 
% Maolin Tian, 2018

global valueTable
curvatureField = zeros(samples.Count, 2);
for s1 = 1:samples.Count
    n1 = samples.tree_ind(s1);
    weigthSum = 0;
    for s2 = cell2mat(tree.sample_ind(tree.ngbr{n1}))'
        n2 = samples.tree_ind(s2);
        d = tree.depth(n2);
        dp = samples.Location(s2,:) - samples.Location(s1,:);
        dx = 1+round((dp(1) - valueTable{d}(1,1)) * 2^(d+2));
        dy = 1+round((dp(2) - valueTable{d}(1,1)) * 2^(d+2));
        Len = size(valueTable{d},1);
        if dx <= 0 || dx > Len || dy <= 0 || dy > Len
            continue;
        end
        F_i = valueTable{d}(dx,2) * valueTable{d}(dy,2) / (tree.width(n2))^2;
        weigthSum = weigthSum + F_i / locationWeights(s2);
        F_i = F_i / locationWeights(s2) * samples.Normal(s2, :);
        curvatureField(s1, :) = curvatureField(s1, :) + F_i;
    end
    curvatureField(s1, 1) = norm(curvatureField(s1, :)) / weigthSum;
end
curvatureField = curvatureField(:, 1);
