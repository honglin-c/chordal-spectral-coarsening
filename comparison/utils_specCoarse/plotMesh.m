function plotFunction(V,F,f, drawEdge)

if nargin < 3
    f = zeros(size(V,1),1);
    EdgeColorBool = false;
elseif nargin == 3
    EdgeColorBool = false;
elseif nargin == 4
    EdgeColorBool = drawEdge;
end
if size(f,1) ~= size(V,1)
    f = zeros(size(V,1),1);
end
    
CMap = myColorMap();
% CMap = gray;
t = tsurf(F,V);
backColor = [1,1,1];
axis equal
axis off

camlight;
if EdgeColorBool == false
    t.EdgeColor = 'none';
else 
    t.EdgeColor = 'black';
end
set(t,fphong, 'FaceVertexCData', f);
set(t, fsoft);
set(gca, 'Visible', 'off')
set(gcf, 'Color', backColor)
% addToolbarExplorationButtons(gcf)
colormap(CMap)
