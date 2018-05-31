function [error, dist] = getError(P1, P2)
%getError Estimate error from P2 to P1(vertex of a isosurface)
% P1 and P2 are N*3 points set. 
% Because dist is the distence from P2 to surface represented by P1, P1 should
% be continues vertex in a isosurface.
P1 = unique(P1, 'rows');

kdObj = KDTreeSearcher(P1);
[match, dist] = knnsearch(kdObj, P2, 'K', 3);
dist = dist(:,1);
% dot(a,corss(b,c))
N = cross((P1(match(:,1),:) - P1(match(:,2),:))', ...
    (P1(match(:,1),:) - P1(match(:,3),:))' );
N = N ./ sqrt(N(:,1).^2 + N(:,2).^2 + N(:,3).^2);
L = P2 - P1(match(:,1),:);
L = (abs(dot(L', N)))';
dist = min(L, dist);
error = sqrt(mean(dist.^2));
disp('error =')
disp(error)
