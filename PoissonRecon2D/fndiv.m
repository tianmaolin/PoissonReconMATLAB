function pp = fndiv(pp, divisor)
%fndiv divide pp.breaks
% suppose pp.order <= 3 and breaks add 1 one time

% suppose pp.order = 3;
if pp.order < 3
    pp.order = 3;
    pp.coefs = [zeros(pp.pieces, 1), pp.coefs];
end

pp.pieces = pp.pieces * divisor;
pp.breaks = linspace(pp.breaks(1), pp.breaks(end), pp.pieces + 1);
new_coefs = zeros(pp.pieces, pp.order);
for k = 1:pp.pieces
    if mod(k - 1, divisor) == 0
        new_coefs(k,:) = pp.coefs(ceil(k/divisor), :);
        n = ceil(k/divisor);
    else
        t = mod(k - 1, divisor) / divisor;
        new_coefs(k,1) = pp.coefs(n, 1);
        new_coefs(k,2) = 2 * pp.coefs(n, 1) * t + pp.coefs(n, 2);
        new_coefs(k,3) = pp.coefs(n, 1) * t * t + pp.coefs(n, 2) * t + pp.coefs(n, 3);
    end
end
pp.coefs = new_coefs;
