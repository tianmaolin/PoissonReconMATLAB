function [tree, samples] = setTree(samples, maxDepth, minDepth)

% TODO: use classdef to save memory
location = samples.Location;
N = 2^maxDepth;
w = 1/N;
u = ceil(double(location(:,1) / w));
v = ceil(double(location(:,2) / w));
A = sparse(u,v,1);
A(A>1) = 1;
A = full(A);
A(N,N)=0;
S = qtdecomp(A,0.5,[1,2^(maxDepth - minDepth)]);

% blocks = repmat(uint8(0),size(S));
% for dim = [512 256 128 64 32 16 8 4 2 1];    
%   numblocks = length(find(S==dim));    
%   if (numblocks > 0)        
%     values = repmat(uint8(1),[dim dim numblocks]);
%     values(2:dim,2:dim,:) = 0;
%     blocks = qtsetblk(blocks,S,dim,values);
%   end
% end
% 
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
    tree.center(:,1)-tree.width/2 == 0 | tree.center(:,2)-tree.width/2 == 0) = true;

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
    ngbr = find(tree.center(k,1)-tree.center(:,1) > -1.5*curW ...
        & tree.center(k,1)-tree.center(:,1) <= 1.5*curW ...
        & tree.center(k,2)-tree.center(:,2) > -1.5*curW ...
        & tree.center(k,2)-tree.center(:,2) <= 1.5*curW);
    tree.ngbr{k} = ngbr;
    tree.ngbr_sample{k} = tree.sample_ind{ngbr};
end

end