function plotSpecParams(moments,momentsSpec,figdir,project,showPlot)

lslopeOrig=momentsSpec.lslope;
rslopeOrig=momentsSpec.rslope;
momentsSpec.lslope(:,moments.elevation>0)=-rslopeOrig(:,moments.elevation>0);
momentsSpec.rslope(:,moments.elevation>0)=-lslopeOrig(:,moments.elevation>0);

levelOrig=momentsSpec.level;
revelOrig=momentsSpec.revel;
momentsSpec.level(:,moments.elevation>0)=-revelOrig(:,moments.elevation>0);
momentsSpec.revel(:,moments.elevation>0)=-levelOrig(:,moments.elevation>0);

lpvelOrig=momentsSpec.lpvel;
rpvelOrig=momentsSpec.rpvel;
momentsSpec.lpvel(:,moments.elevation>0)=-rpvelOrig(:,moments.elevation>0);
momentsSpec.rpvel(:,moments.elevation>0)=-lpvelOrig(:,moments.elevation>0);

momentsSpec.lpvel(isnan(momentsSpec.rpvel))=nan;
momentsSpec.rpvel(isnan(momentsSpec.lpvel))=nan;

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

climsLrwidth=[0,13];
climsLslope=[0,20];
climsRslope=[-20,0];
climsVel=[-15,15];

col1=cat(1,[0,0,0],jet);
col1r=flipud(col1);
col2=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1400 1200],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsSpec.lrwidth(isnan(momentsSpec.lrwidth))=-99;
momentsSpec.lrwidth(isinf(momentsSpec.lrwidth))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lrwidth,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLrwidth);
s1.Colormap=col1;
colorbar
grid on
box on
title('Left edge to right edge width (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpec.lslope(isnan(momentsSpec.lslope))=-99;
momentsSpec.lslope(isinf(momentsSpec.lslope))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLslope);
s3.Colormap=col1;
colorbar
grid on
box on
title('Left slope (dB s m^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

momentsSpec.rslope(isnan(momentsSpec.rslope))=99;
momentsSpec.rslope(isinf(momentsSpec.rslope))=99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.rslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsRslope);
s4.Colormap=col1r;
colorbar
grid on
box on
title('Right slope (dB s m^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s5=nexttile(5);

momentsSpec.level(isnan(momentsSpec.level))=-99;
momentsSpec.level(isinf(momentsSpec.level))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.level,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s5.Colormap=col2;
colorbar
grid on
box on
title('Left edge velocity (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s6=nexttile(6);

momentsSpec.revel(isnan(momentsSpec.revel))=-99;
momentsSpec.revel(isinf(momentsSpec.revel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.revel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s6.Colormap=col2;
colorbar
grid on
box on
title('Right edge velocity (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s7=nexttile(7);

momentsSpec.lpvel(isnan(momentsSpec.lpvel))=-99;
momentsSpec.lpvel(isinf(momentsSpec.lpvel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s7.Colormap=col2;
colorbar
grid on
box on
title('Left peak velocity (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s8=nexttile(8);

momentsSpec.rpvel(isnan(momentsSpec.rpvel))=-99;
momentsSpec.rpvel(isinf(momentsSpec.rpvel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.rpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s8.Colormap=col2;
colorbar
grid on
box on
title('Right peak velocity (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_specParams_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end