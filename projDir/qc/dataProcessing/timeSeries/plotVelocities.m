function plotVelocities(velDual,moments,figdir,project,showPlot)
velBase=velDual(:,:,1);
velHigh=velDual(:,:,2);
velLow=velDual(:,:,3);

velHighFilled=velHigh;
velHighFilled(isnan(velHigh))=velBase(isnan(velHigh));
velLowFilled=velLow;
velLowFilled(isnan(velLow))=velBase(isnan(velLow));

aslGood=moments.asl(~isnan(velBase))./1000;
ylims=[0,max(aslGood)+0.5];

f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(3,2,1);

moments.vel(isnan(moments.vel))=-99;

colTwo=cat(1,[0,0,0],velCols);

hold on
surf(moments.time,moments.asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8 8]);
s1.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity time domain (m s^{-1})')

s2=subplot(3,2,2);

velBase(isnan(velBase))=-99;

hold on
surf(moments.time,moments.asl./1000,velBase,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8 8]);
s2.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity base (m s^{-1})')

s3=subplot(3,2,3);

velHigh(isnan(velHigh))=-99;

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

s4=subplot(3,2,4);

velLow(isnan(velLow))=-99;

hold on
surf(moments.time,moments.asl./1000,velLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8.001 8]);
s4.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity low (m s^{-1})')

s5=subplot(3,2,5);

velHighFilled(isnan(velHighFilled))=-99;

hold on
surf(moments.time,moments.asl./1000,velHighFilled,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8.001 8]);
s5.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity high filled (m s^{-1})')

s6=subplot(3,2,6);

velLowFilled(isnan(velLowFilled))=-99;

hold on
surf(moments.time,moments.asl./1000,velLowFilled,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-8.001 8]);
s6.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity low filled (m s^{-1})')

linkaxes([s1 s2 s3 s4 s5 s6],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_vels100hz_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end