function plotMoments(outString,moments,showPlot,timeBeams,asl,ylimits,figdir,project)

f1 = figure('Position',[200 500 1300 1200],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

s1=subplot(3,2,1);

hold on
surf(timeBeams,asl./1000,moments.powerDB,'edgecolor','none');
view(2);
ylabel('km');
caxis([-110 -70]);
ylim(ylimits);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Power (dB)')

s2=subplot(3,2,2);

hold on
surf(timeBeams,asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('km');
caxis([-40 30]);
ylim(ylimits);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Reflectivity (dBZ)')

s3=subplot(3,2,4);

hold on
surf(timeBeams,asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('km');
caxis([-16 16]);
ylim(ylimits);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity (m s^{-1})')

s4=subplot(3,2,5);

hold on
surf(timeBeams,asl./1000,moments.width,'edgecolor','none');
view(2);
ylabel('km');
caxis([0 2]);
ylim(ylimits);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Spectrum width (m s^{-1})')

s5=subplot(3,2,3);

hold on
surf(timeBeams,asl./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('km');
caxis([-20 50]);
ylim(ylimits);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('SNR (dB)')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_',outString,'_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end