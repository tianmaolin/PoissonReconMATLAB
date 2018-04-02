function ptCloud2D = pointCloud2D(location, normal)
if size(location) ~= size(normal)
    disp('Size of location is not equal to normal!')
    ptCloud2D.Location = [];
    return;
end
if size(location, 2) ~= 2
    disp('Length of row must be 2 !')
    ptCloud2D.Location = [];
    return;
end
ptCloud2D.Location = location;
ptCloud2D.color = [];
ptCloud2D.Normal = normal;
ptCloud2D.Intensity = [];
ptCloud2D.Count = length(location(:,1));
ptCloud2D.XLimits = [min(location(:,1)), max(location(:,1))];
ptCloud2D.YLimits = [min(location(:,2)), max(location(:,2))];
ptCloud2D.dim = 2;