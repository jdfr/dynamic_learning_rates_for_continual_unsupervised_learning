function [Handle]=Figure2pdf(FileName,SizeX,SizeY)

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperSize',[SizeX SizeY]);
set(gcf,'PaperPosition',[0 0 SizeX SizeY]);
saveas(gcf,sprintf('%s.fig',FileName),'fig');
saveas(gcf,sprintf('%s.pdf',FileName),'pdf');
Handle=gcf;
