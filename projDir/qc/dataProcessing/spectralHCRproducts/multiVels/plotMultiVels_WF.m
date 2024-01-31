function plotMultiVels_WF(momVel,momAsl,momTime,shoulderLow,shoulderHigh,velLayers,figdir,project,showPlot)

ylims=[0,12];

% lmin=min(shoulderLow(:),[],'omitmissing');
% cmin=min([-12,lmin]);
% clims=[cmin-0.001,abs(cmin)];
clims=[-13,13];

f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

momVel(isnan(momVel))=-99;

colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

lowLayer=velLayers(:,:,1);
highLayer=velLayers(:,:,2);

dualPartDiff=highLayer-lowLayer;
shoulderDiff=shoulderHigh-shoulderLow;

% difflim=prctile(abs(dualPartDiff(:)),99.5);
% climsDiff=[0,difflim];
% difflimS=prctile(abs(shoulderDiff(:)),99.5);
% climsDiffS=[0,difflimS];

climsDiffS=[0,9];
climsDiff=[0,5];

hold on
surf(momTime,momAsl./1000,momVel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Velocity time domain (m s^{-1})')

s3=nexttile(3);

shoulderHigh(isnan(shoulderHigh))=-99;

hold on
surf(momTime,momAsl./1000,shoulderHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s3.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Velocity high (m s^{-1})')

s4=nexttile(4);

highLayer(isnan(highLayer))=-99;

hold on
surf(momTime,momAsl./1000,highLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles high (m s^{-1})')

s5=nexttile(5);

shoulderLow(isnan(shoulderLow))=-99;

hold on
surf(momTime,momAsl./1000,shoulderLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s5.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Velocity low (m s^{-1})')

s6=nexttile(6);

lowLayer(isnan(lowLayer))=-99;

hold on
surf(momTime,momAsl./1000,lowLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s6.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles low (m s^{-1})')

s7=nexttile(7);

shoulderDiff(isnan(shoulderDiff))=-99;

hold on
surf(momTime,momAsl./1000,shoulderDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiffS);
s7.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Velocity high-low (m s^{-1})')

s8=nexttile(8);

dualPartDiff(isnan(dualPartDiff))=-99;

hold on
surf(momTime,momAsl./1000,dualPartDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s8.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles high-low (m s^{-1})')

linkaxes([s1 s3 s4 s5 s6 s7 s8],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiVel_',datestr(momTime(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momTime(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end