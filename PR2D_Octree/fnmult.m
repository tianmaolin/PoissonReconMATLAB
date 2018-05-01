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
%
% Maolin Tian, Tongji University, 2018

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
    pp2.breaks = [pp2.breaks pp2.breaks(end)+d_break*(1:d_pieces)];
    pp2.coefs = [pp2.coefs; zeros(d_pieces, pp2.order)];
    pp2.pieces = pp2.pieces - d_pieces;
% elseif d_pieces > 0
%     disp('pp2.pieces > pp1.pieces !')
%     return
end

mpp.form = 'pp';
mpp.breaks = pp1.breaks(b1:end);
for i = b1:length(pp1.breaks) - 1
    mpp.coefs(i-b1+1,:) = conv(pp1.coefs(i,:), pp2.coefs(i-b1+1,:));
end
mpp.pieces = length(pp1.breaks) - b1;
mpp.order = pp1.order + pp2.order - 1;
mpp.dim = 1;
