function plotMoments(outString,moments,showPlot,timeBeams,range,ylimUpper,figdir,project)

f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

s1=subplot(5,1,1);

hold on
surf(timeBeams,range./1000,moments.powerDB,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-110 -70]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Power (dB)')
s1pos=s1.Position;
s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];

s2=subplot(5,1,2);

hold on
surf(timeBeams,range./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-40 30]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Reflectivity (dBZ)')
s2pos=s2.Position;
s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

s3=subplot(5,1,3);

hold on
surf(timeBeams,range./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-5 5]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity (m s^{-1})')
s3pos=s3.Position;
s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

s4=subplot(5,1,4);

hold on
surf(timeBeams,range./1000,moments.width,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([0 2]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Spectrum width (m s^{-1})')
s4pos=s4.Position;
s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];

s5=subplot(5,1,5);

hold on
surf(timeBeams,range./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-20 50]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('SNR (dB)')
s5pos=s5.Position;
s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_',outString,'_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end