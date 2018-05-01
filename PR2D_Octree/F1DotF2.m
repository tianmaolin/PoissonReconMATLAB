function Fmap = F1DotF2(F1, F2, w, d)
%F1DotF2 Compute F1_w * F2_w2(t), w2 = 2^d w.
% t = breaks(1):w1:breaks(end).
% F1(x) = F1(-x) and F2(x) = F2(-x)

F1 = fndiv(F1, 2);
F1 = fnscale(F1, w);
F2 = fndiv(F2, 2^(d+1));
F2 = fnscale(F2, 2^d*w);

maxT = F2.breaks(end) - F1.breaks(1);
T = -maxT:w/2:maxT;
V = zeros(length(T),1);
for i = 1:length(T)
    F2_t = fntrans(F2, T(i));
    F = fnmult(F1, F2_t);
    if F.pieces == 0
        V(i) = 0;
    else
        V(i) = fnval(fnint(F),F.breaks(end));
    end
end
Fmap = [T' V];
