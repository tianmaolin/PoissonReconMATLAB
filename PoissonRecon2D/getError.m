function [error, dist] = getError(P1, P2)
%getError Estimate error from P2 to P1(vertex of a contour)
% P1 and P2 are N*2 points set. 
% Because dist is the distence from P2 to line represented by P1, P1 should
% be continues vertex in a contour.
P1 = unique(P1, 'rows');

kdObj = KDTreeSearcher(P1);
[match, dist] = knnsearch(kdObj, P2, 'K', 2);
dist = dist(:,1);
N = P1(match(:,1),:) - P1(match(:,2),:);
N = [N(:,2), - N(:,1)];
N = N ./ sqrt(N(:,1).^2 + N(:,2).^2);
L = P2 - P1(match(:,1),:);
L = L(:,1) .* N(:,1) + L(:,2) .* N(:,2);
L = abs(L);
dist = min(L, dist);
error = sqrt(mean(dist.^2));
disp('error =')
disp(error)
