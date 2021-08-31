function [tree, samples] = setTree(samples, minDepth, maxDepth, feature)
%setTree Set Tree. Support higher depth on feature points.
%
% Maolin Tian, 2018

if minDepth < 0
  minDepth = 0;
  warning('minDepth < 0 !')
end

% TODO: use classdef to save memory
location = samples.Location;
N = 2^maxDepth;
dD = maxDepth - minDepth;
w = 1 / N;
u = ceil(double(location(:, 1) / w));
v = ceil(double(location(:, 2) / w));
% A = sparse(u, v, 1);
% make tree nodes more around points
if nargin < 4
  off = 1;
  u = [u, u, u + off, u + off, u - off, u, u - off];
  v = [v, v + off, v, v + off, v, v - off, v - off];
  v(u <= 0) = [];
  u(u <= 0) = [];
  u(v <= 0) = [];
  v(v <= 0) = [];
  v(u > N) = [];
  u(u > N) = [];
  u(v > N) = [];
  v(v > N) = [];
  A = sparse(u, v, 1);
  A(A > 0) = 2^(dD - 1);
  A = full(A);
else
  off = 2;
  u = [u, u, u + off, u + off, u - off, u, u - off];
  v = [v, v + off, v, v + off, v, v - off, v - off];
  v(u <= 0) = [];
  u(u <= 0) = [];
  u(v <= 0) = [];
  v(v <= 0) = [];
  v(u > N) = [];
  u(u > N) = [];
  u(v > N) = [];
  v(v > N) = [];
  A = sparse(u, v, 2);
  A(A >= 2) = 2^(dD - 2);
  A = full(A);
  for i = 1:size(feature, 1)
    % Suppose old samples used to get weight is new samples(1:N)
    su = ceil(double(feature(i, 1) / w));
    sv = ceil(double(feature(i, 2) / w));
    A(su, sv) = 2^(dD - 1);
    off = 1;
    A(su+off(off <= N - su), sv) = 2^(dD - 1);
    A(su, sv+off(off <= N - sv)) = 2^(dD - 1);
    A(su-off(off < su), sv) = 2^(dD - 1);
    A(su, sv-off(off < sv)) = 2^(dD - 1);
    %         off = 4;
    %         A(su+off(off<=N-su & off<=N-sv), sv+off(off<=N-su & off<=N-sv)) = 2^(dD - 2);
    %         A(su+off(off<=N-su & off <  sv), sv-off(off<=N-su & off <  sv)) = 2^(dD - 2);
    %         A(su-off(off <  su & off<=N-sv), sv+off(off <  su & off<=N-sv)) = 2^(dD - 2);
    %         A(su-off(off <  su & off <  sv), sv-off(off <  su & off <  sv)) = 2^(dD - 2);
  end
  %     for s = depth3Id
  %         su = ceil(double(location(s,1) / w));
  %         sv = ceil(double(location(s,2) / w));
  %         A(su, sv) = 2^(dD - 1);
  %         off = 1;
  %         A(su+off(off<=N-su), sv) = 2^(dD - 1);
  %         A(su, sv+off(off<=N-sv)) = 2^(dD - 1);
  %         A(su-off(off<su), sv) = 2^(dD - 1);
  %         A(su, sv-off(off<sv)) = 2^(dD - 1);
  % %         off = 1;
  % %         A(su+off(off<=N-su & off<=N-sv), sv+off(off<=N-su & off<=N-sv)) = 2^(dD - 1);
  % %         A(su+off(off<=N-su & off <  sv), sv-off(off<=N-su & off <  sv)) = 2^(dD - 1);
  % %         A(su-off(off <  su & off<=N-sv), sv+off(off <  su & off<=N-sv)) = 2^(dD - 1);
  % %         A(su-off(off <  su & off <  sv), sv-off(off <  su & off <  sv)) = 2^(dD - 1);
  %     end
  %     A(A < 1.71 & A >= 1) = 2^(dD - 1);
  %     A(A < 1.93 & A >= 1.71) = 2^(dD - 2);
  %     A(A <= 2 & A >= 1.93) = 2^(dD - 3);
end
A(N, N) = 0;

% If max(block in A) = 2^(n-1), it split at minDepth+n.
% norm([1,0] + [0,1])/2 = 0.7071, -- pi/2
% norm([1,0] + [sqrt(2)/2, sqrt(2)/2])/2 = 0.9239 -- 3/4*pi

% quadtree decomposition
S = qtdecompModified(A, [], [1, 2^(maxDepth - minDepth)]);
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
ind = find(S > 0);
tree.Count = length(ind);
tree.width = full(S(ind)) * w;
tree.depth = -round(log2(tree.width));
tree.center = [rem(ind - 1, N), ceil(ind / N) - 1] * w + 0.5 * tree.width;
tree.isbound = false(tree.Count, 1);
tree.isbound(tree.center(:, 1)+tree.width/2 == 1 | tree.center(:, 2)+tree.width/2 == 1 | ...
  tree.center(:, 1)-tree.width/2 == 0 | tree.center(:, 2)-tree.width/2 == 0 | ...
  tree.center(:, 1)+tree.width*3/2 == 1 | tree.center(:, 2)+tree.width*3/2 == 1 | ...
  tree.center(:, 1)-tree.width*3/2 == 0 | tree.center(:, 2)-tree.width*3/2 == 0) = true;
% tree.isbound(tree.center(:,1)+tree.width/2 == 1 | tree.center(:,2)+tree.width/2 == 1 | ...
%     tree.center(:,1)-tree.width/2 == 0 | tree.center(:,2)-tree.width/2 == 0) = true;

tree.sample_ind = cell(tree.Count, 1);
samples.tree_ind = zeros(samples.Count, 1);
for k = 1:tree.Count
  curW = tree.width(k);
  sam_i = find(tree.center(k, 1)-location(:, 1) > -curW/2 ...
    & tree.center(k, 1)-location(:, 1) <= curW/2 ...
    & tree.center(k, 2)-location(:, 2) > -curW/2 ...
    & tree.center(k, 2)-location(:, 2) <= curW/2);
  samples.tree_ind(sam_i) = k;
  tree.sample_ind{k} = sam_i;
end

tree.ngbr = cell(tree.Count, 1);
tree.ngbr_sample = cell(tree.Count, 1);
for k = 1:tree.Count
  curW = tree.width(k) + tree.width;
  ngbr = find(tree.center(k, 1)-tree.center(:, 1) > -2*curW ...
    & tree.center(k, 1)-tree.center(:, 1) < 2*curW ...
    & tree.center(k, 2)-tree.center(:, 2) > -2*curW ...
    & tree.center(k, 2)-tree.center(:, 2) < 2*curW);
  tree.ngbr{k} = ngbr;
  tree.ngbr_sample{k} = tree.sample_ind{ngbr};
end

end
