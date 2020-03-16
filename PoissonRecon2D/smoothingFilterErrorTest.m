% Plot triangle
n = 500;
im = ones(n,n);
bound = n/10;
minL = ceil(bound);
maxL = floor(n - bound);
for i = minL : maxL
    im(i, minL : n - i) = 0;
end
imshow(im);

% B-spline kernel vs. Gaussian kernel
depth = 7;
b = bspline([-1.5, -0.5, 0.5, 1.5]);
db = fnder(b);
b_w = fnscale(b, 2^(-depth));
kerL = round(3 * n * (2 ^ - depth));
BKernel = zeros(kerL, kerL);
center = round(kerL/2);
for i = 1 : kerL
    for j = 1 : kerL
        x = (center - [i, j]) / n;
        BKernel(i, j)= fnval(b_w, x(1)) * fnval(b_w, x(2));
    end
end
BKernel = BKernel / sum(BKernel(:));
imB = conv2(im, BKernel, 'same');

% figure,hold on
% imshow(imB)
figure,hold on
contour(im, [.5, .5]);
contour(imB, [.5, .5]);



% Gaussian Smoothing
sigma = 2^-depth;
x = - 1.5 * (2 ^ - depth) : 1/n : 1.5 * (2 ^ - depth);
GKernel = gaussmf(x,[sigma 0]);
GKernel = GKernel' * GKernel;
GKernel = GKernel / sum(GKernel(:));
figure, hold on
mesh(GKernel)
mesh(BKernel)

sigma = n*2^-5;
imGaussian = imgaussfilt(im, sigma);
figure, imshow(imGaussian);



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

pp.coefs(:,1) = pp.coefs(:,1)/scale/scale;
pp.coefs(:,2) = pp.coefs(:,2)/scale;
pp.breaks = scale * pp.breaks;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
