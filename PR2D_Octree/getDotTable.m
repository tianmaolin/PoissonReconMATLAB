function [valueTable, dotTable, dotdTable, ddotdTable] = getDotTable(degree, minDepth, maxDepth)
%getDotTable compute b_w, b_w1 * b_w2, b_w1 * b_w2', b_w1' * b_w2'
% x = k 2^(-maxDepth), k is an integer.

if degree == 2
    b = bspline([-1.5, -0.5, 0.5, 1.5]);
elseif degree == 1
	b = bspline([-1, 0, 1]);
end
db = fnder(b);


valueTable = cell(maxDepth - minDepth, 1);
for d = 1 : maxDepth
    b_w = fnscale(b, 2^(-d));
    x = b_w.breaks(1): 2^(-d-2): b_w.breaks(end);
    valueTable{d} = [x', fnval(b_w, x')];
end

dotTable = cell(maxDepth - minDepth, maxDepth - minDepth);
for d1 = minDepth : maxDepth
    for d2 = minDepth : d1
        dotTable{d1,d2} = F1DotF2(b, b,2^(-d1),d1-d2);
    end
end

dotdTable = cell(maxDepth - minDepth, maxDepth - minDepth);
for d1 = minDepth : maxDepth
    % d2 = minDepth : maxDepth because b_w1 * b_w2' is not symmetrical.
    for d2 = minDepth : maxDepth
        % db * 2^d2 is due to (b_w)' = (b'/w)_w
        % - is due to F1DotF2 compute int F1(t) F2(t-x) dt intead of F1*F2
        % = int F1(t) F2(x-t) dt, and F2 = db is odd function.
        if d2 <= d1
            dotdTable{d1,d2} = F1DotF2(b, fncmb(db, 2^d2), 2^(-d1),d1-d2);
            dotdTable{d1,d2}(:,2) = - dotdTable{d1,d2}(:,2);
        else
            dotdTable{d1,d2} = F1DotF2(fncmb(db, 2^d2), b, 2^(-d2),d2-d1);
        end
    end
end

ddotdTable = cell(maxDepth - minDepth, maxDepth - minDepth);
for d1 = minDepth : maxDepth
    for d2 = minDepth : d1
        ddotdTable{d1,d2} = F1DotF2(fncmb(db, 2^d1), fncmb(db, 2^d2), 2^(-d1),d1-d2);
        ddotdTable{d1,d2}(:,2) = - ddotdTable{d1,d2}(:,2);
    end
end

end
