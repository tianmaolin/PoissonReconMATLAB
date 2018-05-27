function [tree, samples] = setTree(samples, minDepth, maxDepth)

if minDepth < 0
    minDepth = 0;
    warning('minDepth < 0 !')
end

% TODO: use classdef to save memory
% octree decomposition
off = 2^-maxDepth;
% location = samples.Location;
location = [samples.Location;...
    samples.Location + off * samples.Normal;...
    samples.Location - off * samples.Normal];
% location = [samples.Location;...
%     samples.Location + off*[1 0 0];...
%     samples.Location - off*[1 0 0];...
%     samples.Location + off*[0 1 0];...
%     samples.Location - off*[0 1 0];...
%     samples.Location + off*[0 0 1];...
%     samples.Location - off*[0 0 1] ];
location(location(:,1) <= 0 | location(:,1) >= 1 |...
    location(:,2) <= 0 | location(:,2) >= 1 |...
    location(:,3) <= 0 | location(:,3) >= 1 , :) = [];
location = [location; 0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1;1 1 1];
OT = OcTreeModified(location , 'binCapacity', 1, 'maxSize', 2^(-minDepth), 'maxDepth', maxDepth + 1);

% RemoveParents' Bins
parentsBin = unique(OT.BinParents);
parentsBin(1) = [];
OT.Points = OT.Points(1:samples.Count,:);
OT.PointBins = OT.PointBins(1:samples.Count,:);
for s = 1:length(OT.Points)
    OT.PointBins(s) = OT.PointBins(s) - length(find(parentsBin<OT.PointBins(s)));
end
% OT.PointBins = OT.PointBins - ...
%     length(find(repmat(parentsBin,length(OT.Points),1)<OT.PointBins));
OT.BinCount = OT.BinCount - length(parentsBin);
OT.BinBoundaries(parentsBin,:) = [];
OT.BinDepths(parentsBin) = [];
OT.Properties.hasParents = false;

% ocTree -> myTree
tree = struct('maxDepth', maxDepth, 'minDepth', minDepth);
tree.Count = OT.BinCount;
tree.width = OT.BinBoundaries(:,4) - OT.BinBoundaries(:,1);
tree.depth = OT.BinDepths';
tree.center = OT.BinBoundaries(:,1:3) + tree.width / 2;
tree.isbound = false(tree.Count,1);
tree.isbound(tree.center(:,1)+tree.width/2 == 1 | tree.center(:,1)-tree.width/2 == 0 |...
    tree.center(:,2)+tree.width/2 == 1 | tree.center(:,2)-tree.width/2 == 0 |...
    tree.center(:,3)+tree.width/2 == 1 | tree.center(:,3)-tree.width/2 == 0) = true;
samples.tree_ind = OT.PointBins;

% Sort myTree, it may be no useful.
[tree.center, index] = sortrows(tree.center);
tree.width = tree.width(index);
tree.depth = tree.depth(index);
tree.isbound = tree.isbound(index);
for s = 1:samples.Count
    samples.tree_ind(s) = find(samples.tree_ind(s)==index, 1);
end

tree.sample_ind = cell(tree.Count,1);
for k = 1:tree.Count
    tree.sample_ind{k} = find(samples.tree_ind == k);
end

tree.ngbr = cell(tree.Count,1);
% tree.ngbr_sample = cell(tree.Count,1);
for k = 1:tree.Count
    curW = tree.width(k)+tree.width;
    % 2*curW is more accurate because points may not on tree's center in setConstantTerms.m.
    tree.ngbr{k} = find(tree.center(k,1)-tree.center(:,1) > -2*curW ...
        & tree.center(k,1)-tree.center(:,1) < 2*curW ...
        & tree.center(k,2)-tree.center(:,2) > -2*curW ...
        & tree.center(k,2)-tree.center(:,2) < 2*curW ...
        & tree.center(k,3)-tree.center(:,3) > -2*curW ...
        & tree.center(k,3)-tree.center(:,3) < 2*curW);
%     tree.ngbr_sample{k} = tree.sample_ind{ngbr};
end

end
