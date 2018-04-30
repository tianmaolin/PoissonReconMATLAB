% Example 2: Decompose 200 random points into bins of 10 points or
%less,shrunk to minimallly encompass their points, then display.
pts = rand(200,3);
OT = OcTree(pts,'binCapacity',10,'style','weighted');
OT.shrink
figure
boxH = OT.plot;
cols = lines(OT.BinCount);
doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
for i = 1:OT.BinCount
    set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
    doplot3(pts(OT.PointBins==i,:),'.','Color',cols(i,:))
end
axis image, view(3)
