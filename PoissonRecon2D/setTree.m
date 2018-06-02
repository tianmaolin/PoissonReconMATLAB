function [tree, samples] = setTree(samples, normalWeight, minDepth, maxDepth)

% TODO: use classdef to save memory
location = samples.Location;
N = 2^maxDepth;
w = 1/N;
u = ceil(double(location(:,1) / w));
v = ceil(double(location(:,2) / w));
% A = sparse(u, v, 1);
% make tree nodes more around points
off = [1,2,4];
% u = [u,       u, u + off, u + off, u - off,       u, u - off];
% v = [v, v + off,       v, v + off,       v, v - off, v - off];
u = [u, u + 0*off, u +   off, u -   off, u - 0*off];
v = [v, v +   off, v + 0*off, v - 0*off, v -   off];
v(u<=0) = [];u(u<=0) = [];
u(v<=0) = [];v(v<=0) = [];
v(u>N) = [];u(u>N) = [];
u(v>N) = [];v(v>N) = [];
A = sparse(u, v, 2);
A(A>2) = 2;
A = full(A);
for s = 1:samples.Count
    su = ceil(double(location(s,1) / w));
    sv = ceil(double(location(s,2) / w));
    A(su, sv) = min(A(su, sv), 1 + normalWeight(s));
    off = [1,2,4];
    A(su+off(off<=N-su), sv) = min(A(su, sv+off(off<=N-sv)), 1 + normalWeight(s));
    A(su, sv+off(off<=N-sv)) = min(A(su, sv+off(off<=N-sv)), 1 + normalWeight(s));
    A(su-off(off<su), sv) = min(A(su-off(off<su), sv), 1 + normalWeight(s));
    A(su, sv-off(off<sv)) = min(A(su, sv-off(off<sv)), 1 + normalWeight(s));
end
% If max(block in A) = 2^(n-1), it split at minDepth+n.
% norm([1,0] + [0,1])/2 = 0.7071, -- pi/2
% norm([1,0] + [sqrt(2)/2, sqrt(2)/2])/2 = 0.9239 -- 3/4*pi
A(A < 1.71 & A >= 1) = 2^(maxDepth - minDepth - 1);
A(A < 1.93 & A >= 1.71) = 2^(maxDepth - minDepth - 2);
A(A <= 2 & A >= 1.93) = 2^(maxDepth - minDepth - 3);
A(N,N)=0;


% quadtree decomposition
S = qtdecompModified(A, [], [1 2^(maxDepth - minDepth)]);
% S = qtdecomp(A,0.5,[1,2^(maxDepth - minDepth)]);

% blocks = repmat(uint8(0),size(S));
% for dim = [512 256 128 64 32 16 8 4 2 1]    
%   numblocks = length(find(S==dim));    
%   if (numblocks > 0)        
%     values = repmat(uint8(1),[dim dim numblocks]);
%     values(2:dim,2:dim,:) = 0;
%     blocks = qtsetblk(blocks,S,dim,values);
%   end
% end
% blocks(end,1:end) = 1;
% blocks(1:end,end) = 1;
% figure
% imshow(blocks,[])

tree = struct('maxDepth', maxDepth, 'minDepth', minDepth);
ind = find(S>0);
tree.Count = length(ind);
tree.width = full(S(ind)) * w;
tree.depth = -round(log2(tree.width));
tree.center = [rem(ind-1,N),ceil(ind/N)-1] * w + 0.5*tree.width;
tree.isbound = false(tree.Count,1);
tree.isbound(tree.center(:,1)+tree.width/2 == 1 | tree.center(:,2)+tree.width/2 == 1 | ...
    tree.center(:,1)-tree.width/2 == 0 | tree.center(:,2)-tree.width/2 == 0 | ...
    tree.center(:,1)+tree.width*3/2 == 1 | tree.center(:,2)+tree.width*3/2 == 1 | ...
    tree.center(:,1)-tree.width*3/2 == 0 | tree.center(:,2)-tree.width*3/2 == 0) = true;
% tree.isbound(tree.center(:,1)+tree.width/2 == 1 | tree.center(:,2)+tree.width/2 == 1 | ...
%     tree.center(:,1)-tree.width/2 == 0 | tree.center(:,2)-tree.width/2 == 0) = true;

tree.sample_ind = cell(tree.Count,1);
samples.tree_ind = zeros(samples.Count,1);
for k = 1:tree.Count
    curW = tree.width(k);
    sam_i = find(tree.center(k,1)-location(:,1) > -curW/2 ...
        & tree.center(k,1)-location(:,1) <= curW/2 ...
        & tree.center(k,2)-location(:,2) > -curW/2 ...
        & tree.center(k,2)-location(:,2) <= curW/2);
    samples.tree_ind(sam_i) = k;
    tree.sample_ind{k} = sam_i;
end

tree.ngbr = cell(tree.Count,1);
tree.ngbr_sample = cell(tree.Count,1);
for k = 1:tree.Count
    curW = tree.width(k)+tree.width;
    ngbr = find(tree.center(k,1)-tree.center(:,1) > -2*curW ...
        & tree.center(k,1)-tree.center(:,1) < 2*curW ...
        & tree.center(k,2)-tree.center(:,2) > -2*curW ...
        & tree.center(k,2)-tree.center(:,2) < 2*curW);
    tree.ngbr{k} = ngbr;
    tree.ngbr_sample{k} = tree.sample_ind{ngbr};
end

end
