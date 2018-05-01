function ptCloudNormalized = pcdownsample2D(ptCloud, width)
%pcdownsample2D Downsample the 2-D ptCloud
location = [ptCloud.Location,zeros(ptCloud.Count,1)];
normal = [ptCloud.Normal,zeros(ptCloud.Count,1)];

p = pointCloud(location,'Normal',normal);
p = pcdownsample(p,'gridAverage',width);

location = p.Location(:,1:2);
normal = p.Normal(:,1:2);

ptCloudNormalized = pointCloud2D(location, normal);
end
