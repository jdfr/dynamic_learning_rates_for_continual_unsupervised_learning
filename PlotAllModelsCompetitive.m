function PlotAllModelsCompetitive(do_comp);

ColorSpaces = {'HSL','HSV','Lab','Luv','RGB','YCbCr'};
NumColorSpaces = length(ColorSpaces);
for NdxColorSpace=1:NumColorSpaces,
  MyColorSpace = ColorSpaces{NdxColorSpace};
  load(['./Modelos_' MyColorSpace '.mat'],'Modelos');
  if do_comp
    h=PlotCompetitive(Modelos);
    savefig([MyColorSpace '_competitive.fig']);
    close(h);
    h=PlotCompetitivePCA(Modelos);
    savefig([MyColorSpace '_competitive_pca.fig']);
    close(h);
  else
    h=PlotSOM(Modelos);
    savefig([MyColorSpace '_SOM.fig']);
    close(h);
    h=PlotSOMPCA(Modelos);
    savefig([MyColorSpace '_SOM_pca.fig']);
    close(h);
  endif
end

