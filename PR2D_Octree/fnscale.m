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

pp.coefs(:,1) = pp.coefs(:,1)/scale/scale;
pp.coefs(:,2) = pp.coefs(:,2)/scale;
pp.breaks = scale * pp.breaks;
