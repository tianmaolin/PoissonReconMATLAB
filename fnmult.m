function mpp = fnmult(pp1, pp2)
%fnmult Multiply function
% mpp = pp1 * pp2. pp1 defined on [b3,b4] and pp2 defined on [b5,b6]. b4 -
% b3 = b6 - b5¡£ [b1, b2] = [b3,b4] intersect [b5,b6]. And breaks add 1 one
% time. isint(b3) && isint(b5) = 1.

% v1,v2 in conv(v1, v2) must be a vector, so the following command failed.
% mpp = spmak(fnbrk(pp1,'b'), conv(fnbrk(pp1,'c'), fnbrk(pp2,'c')));

b1 = find(pp1.breaks == pp2.breaks(1), 1);
if isempty(b1)
    pp_temp = pp1;
    pp1 = pp2;
    pp2 = pp_temp;

    b1 = find(pp1.breaks == pp2.breaks(1), 1);
    if isempty(b1)
        mpp.pieces = 0;
        return;
    end
end

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

mpp.form = 'pp';
mpp.breaks = pp1.breaks(b1:end);
for i = b1:length(pp1.breaks) - 1
    mpp.coefs(i-b1+1,:) = conv(pp1.coefs(i,:), pp2.coefs(i-b1+1,:));
end
mpp.pieces = length(pp1.breaks) - b1;
mpp.order = pp1.order + pp2.order - 1;
mpp.dim = 1;
