% function Sum = getLinearBsplneSum(x, N, p_N, B_N)
% 
% if floor(p_N / N) ~= p_N / N
%     Sum = 0;
%     return;
% end
% 
% w = 1/N;
% p = 0:1/p_N:1;
% [px, py] = meshgrid(p, p);
% px = px(:);
% py = py(:);
% map_xy = @(n) [rem(n, N),floor(n / N) + 1];
% 
% neighbor = cell(N^2,1);
% for n = 1:N^2
%     o = (map_xy(n) - 0.5) * w;
%     a = find(abs(px(:) - o(1)) < 2 * w);
%     b = find(abs(py(:) - o(2)) < 2 * w);
%     neighbor{n} = intersect(a,b)';
% end
% 
% % B = bspline([-1.5,-0.5,0.5,1.5]); % basic function
% B = bspline([-1,0,1]); % basic function
% F = zeros(B_N, B_N);
% for i = 1:B_N
%     for j = 1:B_N
%         Bx = (i-1) / B_N * 2;
%         By = (j-1) / B_N * 2;
%         if Bx > B.breaks(end) || By > B.breaks(end)
%             F(i,j) = 0;
%         else
%             F(i,j) = fnval(B, Bx) * fnval(B, By) / w^2;
%         end
%     end
% end
% 
% map_xy = @(n) [rem(n, p_N),floor(n / p_N) + 1];
% Sum = zeros(p_N ^ 2,1);
% for n = 1:p_N^2
%     for i = neighbor{n}
%         c = (map_xy(i) - 0.5) * w;
%         t = [abs(px(map_xy(n)) - c(1)) / w, abs(px(map_xy(n)) - c(2)) / w];
%         t = [floor(t(1) * p_N/2) + 1, floor(t(2) * p_N/2) + 1];
%         if max(t) > p_N
%             continue;
%         end
%         B = F(t(1)) * F(t(2)) / w^2;
%         Sum(n) = Sum(n) + x(n) * B;
%     end
% end
