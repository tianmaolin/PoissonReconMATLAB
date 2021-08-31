function [valueTable, dotTable, dotdTable, ddotdTable] =...
    valueDotTable(minDepth, maxDepth, degree)
%valueDotTable compute value table b_w, and dot product tables,
% b_w1 * b_w2, b_w1 * b_w2', b_w1' * b_w2'.
% For the means of the functions see the technical document.
%
% degree: the degree of B-spline, i.e., degree=1 is linear function.
%
% Maolin Tian, 2018

if nargin < 3
  degree = 2;
end
if degree == 2
  b = bspline([-1.5, -0.5, 0.5, 1.5]);
elseif degree == 1
  b = bspline([-1, 0, 1]);
end
db = fnder(b);


valueTable = cell(maxDepth-minDepth, 1);
for d = 1:maxDepth
  b_w = fnscale(b, 2^(-d));
  x = b_w.breaks(1):2^(-d - 2):b_w.breaks(end);
  valueTable{d} = [x', fnval(b_w, x')];
end

dotTable = cell(maxDepth-minDepth, maxDepth-minDepth);
for d1 = minDepth:maxDepth
  for d2 = minDepth:d1
    dotTable{d1, d2} = F1DotF2(b, b, 2^(-d1), d1-d2);
  end
end

dotdTable = cell(maxDepth-minDepth, maxDepth-minDepth);
for d1 = minDepth:maxDepth
  % d2 = minDepth : maxDepth because b_w1 * b_w2' is not symmetrical.
  for d2 = minDepth:maxDepth
    % db * 2^d2 is due to (b_w)' = (b'/w)_w
    % - is due to F1DotF2 compute int F1(t) F2(t-x) dt intead of F1*F2
    % = int F1(t) F2(x-t) dt, and F2 = db is odd function.
    if d2 <= d1
      dotdTable{d1, d2} = F1DotF2(b, fncmb(db, 2^d2), 2^(-d1), d1-d2);
      dotdTable{d1, d2}(:, 2) = -dotdTable{d1, d2}(:, 2);
    else
      dotdTable{d1, d2} = F1DotF2(fncmb(db, 2^d2), b, 2^(-d2), d2-d1);
    end
  end
end

ddotdTable = cell(maxDepth-minDepth, maxDepth-minDepth);
for d1 = minDepth:maxDepth
  for d2 = minDepth:d1
    ddotdTable{d1, d2} = F1DotF2(fncmb(db, 2^d1), fncmb(db, 2^d2), 2^(-d1), d1-d2);
    ddotdTable{d1, d2}(:, 2) = -ddotdTable{d1, d2}(:, 2);
  end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fmap = F1DotF2(F1, F2, w, d)
%F1DotF2 Compute F1_w * F2_w2(t), w2 = 2^d w.
% t = breaks(1):w1:breaks(end).
% Suppose F1(x) = F1(-x) and F2(x) = F2(-x)

F1 = fndiv(F1, 2);
F1 = fnscale(F1, w);
F2 = fndiv(F2, 2^(d + 1));
F2 = fnscale(F2, 2^d*w);

maxT = F2.breaks(end) - F1.breaks(1);
T = -maxT:w / 2:maxT;
V = zeros(length(T), 1);
for i = 1:length(T)
  F2_t = fntrans(F2, T(i));
  F = fnmult(F1, F2_t);
  if F.pieces == 0
    V(i) = 0;
  else
    V(i) = fnval(fnint(F), F.breaks(end));
  end
end
Fmap = [T', V];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pp = fndiv(pp, divisor)
%fndiv divide pp.breaks
% suppose pp.order <= 3 and breaks add 1 one time

% suppose pp.order = 3;
if pp.order == 2
  pp.order = 3;
  pp.coefs = [zeros(pp.pieces, 1), pp.coefs];
elseif pp.order == 1
  pp.order = 3;
  pp.coefs = [zeros(pp.pieces, 1), zeros(pp.pieces, 1), pp.coefs];
elseif pp.order == 0
  pp.order = 3;
  pp.coefs = zeros(pp.pieces, 3);
end

pp.pieces = pp.pieces * divisor;
pp.breaks = linspace(pp.breaks(1), pp.breaks(end), pp.pieces+1);
new_coefs = zeros(pp.pieces, pp.order);
for k = 1:pp.pieces
  if mod(k-1, divisor) == 0
    new_coefs(k, :) = pp.coefs(ceil(k / divisor), :);
    n = ceil(k/divisor);
  else
    t = mod(k-1, divisor) / divisor;
    new_coefs(k, 1) = pp.coefs(n, 1);
    new_coefs(k, 2) = 2 * pp.coefs(n, 1) * t + pp.coefs(n, 2);
    new_coefs(k, 3) = pp.coefs(n, 1) * t * t + pp.coefs(n, 2) * t + pp.coefs(n, 3);
  end
end
pp.coefs = new_coefs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pp = fnscale(pp, scale)
%fnscale pp(x) = pp(x/scale)
% pp.order should <= 3

% suppose pp.order = 3;
if pp.order == 2
  pp.order = 3;
  pp.coefs = [zeros(pp.pieces, 1), pp.coefs];
elseif pp.order == 1
  pp.order = 3;
  pp.coefs = [zeros(pp.pieces, 1), zeros(pp.pieces, 1), pp.coefs];
elseif pp.order == 0
  pp.order = 3;
  pp.coefs = zeros(pp.pieces, 3);
end

pp.coefs(:, 1) = pp.coefs(:, 1) / scale / scale;
pp.coefs(:, 2) = pp.coefs(:, 2) / scale;
pp.breaks = scale * pp.breaks;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tpp = fntrans(pp1, t)
%fnmult Translate function
tpp = pp1;
tpp.breaks = tpp.breaks + t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mpp = fnmult(pp1, pp2)
%fnmult Multiply function
% d_breaks of pp1 and pp2 should be equal
% v1,v2 in conv(v1, v2) must be a vector, so the following command failed.
% mpp = spmak(fnbrk(pp1,'b'), conv(fnbrk(pp1,'c'), fnbrk(pp2,'c')));
%
% b2 = length(pp2.breaks) - b1 + 1;
%             b1
%             |
%             |
%            \ /
% pp1:    -----------
% pp2:        -----------
%                  / \
%                   |
%                   |
%                   b2


b1 = find(abs(pp1.breaks - pp2.breaks(1)) < 1e-13, 1);
if isempty(b1)
  pp_temp = pp1;
  pp1 = pp2;
  pp2 = pp_temp;

  b1 = find(abs(pp1.breaks - pp2.breaks(1)) < 1e-13, 1);
  if isempty(b1)
    mpp.pieces = 0;
    return;
  end
end

% suppose pp2 has the same length with pp1
d_pieces = pp1.pieces - pp2.pieces;
if d_pieces > 0
  d_break = pp2.breaks(2) - pp2.breaks(1);
  pp2.breaks = [pp2.breaks, pp2.breaks(end) + d_break * (1:d_pieces)];
  pp2.coefs = [pp2.coefs; zeros(d_pieces, pp2.order)];
  pp2.pieces = pp2.pieces - d_pieces;
  % elseif d_pieces > 0
  %     disp('pp2.pieces > pp1.pieces !')
  %     return
end

mpp.form = 'pp';
mpp.breaks = pp1.breaks(b1:end);
for i = b1:length(pp1.breaks) - 1
  mpp.coefs(i-b1+1, :) = conv(pp1.coefs(i, :), pp2.coefs(i - b1 + 1, :));
end
mpp.pieces = length(pp1.breaks) - b1;
mpp.order = pp1.order + pp2.order - 1;
mpp.dim = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
