function plotBasicMoments(moments,cf,plotTimeAll,figdir,project,showPlot)

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

climsPow=[-105,-35];
climsDbz=[-40,30];
climsSnr=[-20,60];
climsNcp=[0,1.2];
colDiff=cat(1,[0,0,0],jet);

%% Figure
f1 = figure('Position',[200 500 1600 830],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

moments.powerV(isnan(cf.VEL_MASKED))=-999;

hold on
surf(moments.time,moments.asl./1000,moments.powerV,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsPow);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain DBMV (dBm)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s2=nexttile(2);

moments.dbz(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDbz);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain DBZ (dB)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

moments.snr(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSnr);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain SNR (dB)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

moments.ncp(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.ncp,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsNcp);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain NCP (dB)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_basicMoments_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end