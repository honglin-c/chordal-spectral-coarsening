function plotFMap(fMap)
imagesc(fMap,[-max(max(abs(fMap))) max(max(abs(fMap)))])
colormap(gca,myColorMap())
axis image;