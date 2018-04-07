function b = getConstantTerms(grid, samples, weight, div)
%getb b_i = int d B_i * F dp

N = grid.size(1);
% map = @(i,j) i + (j-1)*N;
map_ij = @(n) [rem(n, N),ceil(n / N)];
w = grid.width;

dx = 1 / div;
x = [-3:dx:-dx, 0:dx:3];
ppformInt_c = ceil(length(x) / 2);
% TODO: change X,Y
[X, Y] = meshgrid(x,x);
b_i_x = zeros(length(x), length(x));
for a = 1:length(x)
    for b = 1:length(x)
        b_i_x(a, b) = int_dxB_B([X(a,b), Y(a,b)], div) / w ^ 3;
    end
end

ngbrT = repmat(-3:1:3, 7, 1) + repmat((-3*N:N:3*N)', 1, 7);
grid.neighbor = cell(N * N, 1);
for n = 1 : N * N
    ngbr = ngbrT + n;
    for m = (ngbr(:))'
        if m < 1 || m > N * N
            continue
        end
        if isempty(grid.sample_ind{m})
            continue
        end
        grid.neighbor{n} = [grid.neighbor{n}, (grid.sample_ind{m})'];
    end
end

b = zeros(N^2, 1);
for n = 1:N^2
%     if isempty(grid.neighbor{n})
%         continue;
%     end
    for s = grid.neighbor{n}
        p = samples.Location(s,:);
        normal = samples.Normal(s,:);
        o = (map_ij(n) - 0.5) * grid.width;
        if max(abs(p - o)) >= 3 * w
            continue;
        end
        t_x = (p(1)-o(1))/w; 
        t_y = (p(2)-o(2))/w;
        t_m = round(t_y/dx) + ppformInt_c;
        t_n = round(t_x/dx) + ppformInt_c;
        temp = [b_i_x(t_m,t_n); b_i_x(t_n,t_m)];
        t_n = normal * temp;
        b(n) = b(n) + t_n / weight(s);
    end
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b_i_x = int_dxB_B(trans, div)
% B = bspline([-1.5, -0.5, 0.5, 1.5]); % basic function
B = bspline([-1, 0, 1]); % basic function
B = fndiv(B, div);
dB = fnder(B);
b_i_x = fn_int_F_Ft(B, dB, trans(1)) * fn_int_F_Ft(B, B, trans(2));
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
end