function [ samples , pointsFeature ]= refineTreeAdaptiveToCurvature(ptCloud, samp1, maxDepth, feature)
% Maolin Tian, 2018
location = [];
normal = [];
samW = 2^-maxDepth;
for s = 1:size(feature,1)
    id = ptCloud.Location(:,1) < feature(s,1)+samW & ptCloud.Location(:,1) > feature(s,1)-samW &...
        ptCloud.Location(:,2) < feature(s,2)+samW & ptCloud.Location(:,2) > feature(s,2)-samW &...
        ptCloud.Location(:,3) < feature(s,3)+samW & ptCloud.Location(:,3) > feature(s,3)-samW;
    location = [location; ptCloud.Location(id,:)];
    normal = [normal; ptCloud.Normal(id,:)];
    id = samp1.Location(:,1) < feature(s,1)+samW & samp1.Location(:,1) > feature(s,1)-samW &...
        samp1.Location(:,2) < feature(s,2)+samW & samp1.Location(:,2) > feature(s,2)-samW &...
        samp1.Location(:,3) < feature(s,3)+samW & samp1.Location(:,3) > feature(s,3)-samW;
    samp1.Location(id,:) = [];
    samp1.Normal(id,:) = [];
end
if isempty(location)
    pc2 = pointCloud(location, 'Normal', normal);
else
    pc2 = pcdownsampleConst(pointCloud(location, 'Normal', normal), 2^(-maxDepth));
end
samples = struct('Count', size(samp1.Location,1) + size(pc2.Location,1), 'Location', [samp1.Location;pc2.Location],'Normal', [samp1.Normal;pc2.Normal]);
pointsFeature = [false(size(samp1.Location,1),1); true(size(pc2.Location,1),1)];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ptCloudOut = pcdownsampleConst(ptCloudIn, width)

if ptCloudIn.Count == 0
    ptCloudOut = ptCloudIn;
    return
end

id = ceil(ptCloudIn.Location / width);
[~ , ia] = unique(id,'rows');

ptCloudOut = pointCloud(ptCloudIn.Location(ia,:), 'Normal', ptCloudIn.Normal(ia,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%