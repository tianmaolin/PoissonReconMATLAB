% Example 1: Decompose 200 random points into bins of 20 points or
% less,then display each bin with its points in a separate colour.
pts = pd.Location;
OT = OcTree(pts,'maxDepth',6,'binCapacity',1);
figure
boxH = OT.plot;
cols = lines(OT.BinCount);
doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
for i = 1:OT.BinCount
    set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
    doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
end
axis image, view(3)
