function plotVelocities(velDual,moments,figdir,project,showPlot)
velBase=velDual(:,:,1);
velHigh=velDual(:,:,2);
velLow=velDual(:,:,3);

aslGood=moments.asl(~isnan(velBase))./1000;
ylims=[0,max(aslGood)+0.5];

f1 = figure('Position',[200 500 1800 900],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(2,2,1);

hold on
surf(moments.time,moments.asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8 8]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity time domain (m s^{-1})')

s2=subplot(2,2,2);

hold on
surf(moments.time,moments.asl./1000,velBase,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8 8]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity base (m s^{-1})')

s3=subplot(2,2,3);

velHigh(isnan(velHigh))=-99;

colTwo=cat(1,[0,0,0],velCols);

hold on
surf(moments.time,moments.asl./1000,velHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8.001 8]);
s3.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity high (m s^{-1})')

s4=subplot(2,2,4);

velLow(isnan(velLow))=-99;

hold on
surf(moments.time,moments.asl./1000,velLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8.001 8]);
s4.Colormap=colTwo;
ylim(ylims);
colorbar
grid on
box on
title('Velocity low (m s^{-1})')

linkaxes([s1 s2 s3 s4],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_vels100hz_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end