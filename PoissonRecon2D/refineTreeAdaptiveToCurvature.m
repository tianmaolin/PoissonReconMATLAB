function samples = refineTreeAdaptiveToCurvature(ptCloud2d, maxDepth, pc, feature)
%getNormalWeight normalWeights(q) = sum_s F_s(q) / s.locationWeights * s.N
%/ sum_s F_s(q) / s.locationWeights
%
% v2 (original Poisson recon version) depth 7 = v3 (curvature-adaptive 
% Poisson recon version) (p=0 and depth=8) or (p=1 and depth=7)
%
% Maolin Tian, 2018

% feature = pc.Location(normalWeights < 2,:);
% norm([1,0] + [sqrt(2)/2, sqrt(2)/2])/2 = 0.9239 --- 3/4*pi
location = [];
normal = [];
samW = 2^-maxDepth;
for s = 1:size(feature,1)
    id = ptCloud2d.Location(:,1) < feature(s,1)+samW & ptCloud2d.Location(:,1) > feature(s,1)-samW &...
        ptCloud2d.Location(:,2) < feature(s,2)+samW & ptCloud2d.Location(:,2) > feature(s,2)-samW;
    location = [location; ptCloud2d.Location(id,:)];
    normal = [normal; ptCloud2d.Normal(id,:)];
    id = pc.Location(:,1) < feature(s,1)+samW & pc.Location(:,1) > feature(s,1)-samW &...
        pc.Location(:,2) < feature(s,2)+samW & pc.Location(:,2) > feature(s,2)-samW;
    pc.Location(id,:) = [];
    pc.Normal(id,:) = [];
end
pc2 = pcdownsample2D(pointCloud2D(location,normal), 2^-maxDepth);
samples = pointCloud2D([pc.Location;pc2.Location], [pc.Normal;pc2.Normal]);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudNormalized = pcdownsample2D(ptCloud, width)
%pcdownsample2D Downsample the 2-D ptCloud

if ptCloud.Count == 0
    ptCloudNormalized = ptCloud;
    return
end

id = ceil(ptCloud.Location / width);
[~ , ia] = unique(id,'rows');

ptCloudNormalized = pointCloud2D(ptCloud.Location(ia,:), ptCloud.Normal(ia,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%