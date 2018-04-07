function mpp = fnmult(pp1, pp2)
%fnmult Multiply function
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

mpp.form = 'pp';
mpp.breaks = pp1.breaks(b1:end);
for i = b1:length(pp1.breaks) - 1
    mpp.coefs(i-b1+1,:) = conv(pp1.coefs(i,:), pp2.coefs(i-b1+1,:));
end
mpp.pieces = length(pp1.breaks) - b1;
mpp.order = pp1.order + pp2.order - 1;
mpp.dim = 1;
