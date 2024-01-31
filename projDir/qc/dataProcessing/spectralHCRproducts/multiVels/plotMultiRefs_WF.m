function plotMultiRefs_WF(momDbz,momAsl,momTime,shoulderLowDbz,shoulderHighDbz,dbzLayers,figdir,project,showPlot)

ylims=[0,12];

clims=[-20,25];

colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

lowLayer=dbzLayers(:,:,1);
highLayer=dbzLayers(:,:,2);

dualPartDiff=highLayer-lowLayer;
shoulderDiffDbz=shoulderHighDbz-shoulderLowDbz;

% difflim=prctile(abs(dualPartDiff(:)),99.5);
% climsDiff=[-difflim,difflim];
% difflimS=prctile(abs(shoulderDiffDbz(:)),99.5);
% climsDiffS=[-difflimS,difflimS];

climsDiffS=[-20,20];
climsDiff=[-17,17];

f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

momDbz(isnan(momDbz))=-99;

hold on
surf(momTime,momAsl./1000,momDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Reflectivity time domain (dBZ)')

s3=nexttile(3);

shoulderHighDbz(isnan(shoulderHighDbz))=-99;

hold on
surf(momTime,momAsl./1000,shoulderHighDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s3.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Reflectivity high (dBZ)')

s4=nexttile(4);

highLayer(isnan(highLayer))=-99;

hold on
surf(momTime,momAsl./1000,highLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles high (dBZ)')

s5=nexttile(5);

shoulderLowDbz(isnan(shoulderLowDbz))=-99;

hold on
surf(momTime,momAsl./1000,shoulderLowDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s5.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Reflectivity low (dBZ)')

s6=nexttile(6);

lowLayer(isnan(lowLayer))=-99;

hold on
surf(momTime,momAsl./1000,lowLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s6.Colormap=colDiff;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles low (dBZ)')

s7=nexttile(7);

shoulderDiffDbz(isnan(shoulderDiffDbz))=-99;

hold on
surf(momTime,momAsl./1000,shoulderDiffDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiffS);
s7.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Reflectivity high-low (dB)')

s8=nexttile(8);

dualPartDiff(isnan(dualPartDiff))=-99;

hold on
surf(momTime,momAsl./1000,dualPartDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s8.Colormap=colTwo;
ylim(ylims);
xlim([momTime(1),momTime(end)]);
colorbar
grid on
box on
title('Dual particles high-low (dB)')

linkaxes([s1 s3 s4 s5 s6 s7 s8],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiRef_',datestr(momTime(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momTime(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end