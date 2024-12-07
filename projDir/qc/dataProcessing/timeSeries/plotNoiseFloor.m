function plotNoiseFloor(moments,noiseFloor,cf,plotTimeAll,figdir,project,showPlot)

aslGood=moments.asl(~isnan(moments.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsNoiseFloor=[-65,-45];
colDiff=cat(1,[0,0,0],jet);

%% Figure
f1 = figure('Position',[200 500 800 415],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

noiseFloor(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,noiseFloor,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsNoiseFloor);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('NoiseFloor (dB)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);


set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_noiseFloor_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end