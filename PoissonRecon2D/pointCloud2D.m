function ptCloud2D = pointCloud2D(location, normal)
%pointCloud2D Struct for storing a 2-D point cloud.
%
% Maolin Tian, 2018
if size(location) ~= size(normal)
  error('Size of location is not equal to normal!')
end
if isempty(location)
  ptCloud2D = struct('Location', [], 'Normal', [], 'Count', 0, 'XLimits', [], 'YLimits', [], 'dim', 2);
  return
end
if size(location, 2) ~= 2
  error('The number of columns must be 2 !')
end

ptCloud2D.Location = location;
ptCloud2D.Normal = normal;
ptCloud2D.Count = length(location(:, 1));
ptCloud2D.XLimits = [min(location(:, 1)), max(location(:, 1))];
ptCloud2D.YLimits = [min(location(:, 2)), max(location(:, 2))];
ptCloud2D.dim = 2;
end
