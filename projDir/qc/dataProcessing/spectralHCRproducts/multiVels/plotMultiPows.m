function plotMultiPows(moments,dataCF,peakLow,peakHigh,figdir,project,showPlot)
moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

% lmin=min(shoulderLow(:),[],'omitmissing');
% cmin=min([-12,lmin]);
% clims=[cmin-0.001,abs(cmin)];
clims=[-120,-20];
clims2=[-80,20];

f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

moments.powerV(isnan(moments.vel))=-999;

colDiff=cat(1,[0,0,0],jet);

dualPartDiff=peakHigh-peakLow;

% difflim=prctile(abs(dualPartDiff(:)),99.5);
% climsDiff=[0,difflim];
% difflimS=prctile(abs(shoulderDiff(:)),99.5);
% climsDiffS=[0,difflimS];

climsDiffS=[0,12];
climsDiff=[-10,10];

hold on
surf(moments.time,moments.asl./1000,moments.powerV,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Power time domain (dB)')

% s2=nexttile(2);
% 
% dataCF.VEL_MASKED(isnan(dataCF.VEL_MASKED))=-99;
% 
% hold on
% surf(dataCF.time,dataCF.asl./1000,dataCF.VEL_MASKED,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s2.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity time domain cfRadial (m s^{-1})')

% s3=nexttile(3);
% 
% shoulderHigh(isnan(shoulderHigh))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,shoulderHigh,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s3.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity high (m s^{-1})')

s4=nexttile(4);

peakHigh(isnan(peakHigh))=-999;

hold on
surf(moments.time,moments.asl./1000,peakHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims2);
s4.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles high (dB)')

% s5=nexttile(5);
% 
% shoulderLow(isnan(shoulderLow))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,shoulderLow,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s5.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity low (m s^{-1})')

s6=nexttile(6);

peakLow(isnan(peakLow))=-999;

hold on
surf(moments.time,moments.asl./1000,peakLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims2);
s6.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles low (m s^{-1})')

% s7=nexttile(7);
% 
% shoulderDiff(isnan(shoulderDiff))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,shoulderDiff,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(climsDiffS);
% s7.Colormap=colDiff;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity high-low (m s^{-1})')

s8=nexttile(8);

dualPartDiff(isnan(dualPartDiff))=-99;

hold on
surf(moments.time,moments.asl./1000,dualPartDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s8.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles high-low (m s^{-1})')

%linkaxes([s1 s3 s4 s5 s6 s7 s8],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiVelPow_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end