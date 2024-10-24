function plotMomentsSpecParams(momentsSpec,figdir,project,showPlot)

momentsSpec.velRaw(:,momentsSpec.elevation>0)=-momentsSpec.velRaw(:,momentsSpec.elevation>0);
momentsSpec.skew(:,momentsSpec.elevation>0)=-momentsSpec.skew(:,momentsSpec.elevation>0);

lslopeOrig=momentsSpec.lslope;
rslopeOrig=momentsSpec.rslope;
momentsSpec.lslope(:,momentsSpec.elevation>0)=-rslopeOrig(:,momentsSpec.elevation>0);
momentsSpec.rslope(:,momentsSpec.elevation>0)=-lslopeOrig(:,momentsSpec.elevation>0);

levelOrig=momentsSpec.level;
revelOrig=momentsSpec.revel;
momentsSpec.level(:,momentsSpec.elevation>0)=-revelOrig(:,momentsSpec.elevation>0);
momentsSpec.revel(:,momentsSpec.elevation>0)=-levelOrig(:,momentsSpec.elevation>0);

lpvelOrig=momentsSpec.lpvel;
rpvelOrig=momentsSpec.rpvel;
momentsSpec.lpvel(:,momentsSpec.elevation>0)=-rpvelOrig(:,momentsSpec.elevation>0);
momentsSpec.rpvel(:,momentsSpec.elevation>0)=-lpvelOrig(:,momentsSpec.elevation>0);

aslGood=momentsSpec.asl(~isnan(momentsSpec.skew))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,3];
climsSkew=[-3,3];
climsKurt=[-6,6];
climsLrwidth=[0,13];
climsLslope=[0,20];
climsRslope=[-20,0];
climsVel=[-15,15];

col1=cat(1,[0,0,0],jet);
col1r=flipud(col1);
col2=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 2100 1200],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(4,3,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsSpec.width(isnan(momentsSpec.width))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=col1;
colorbar
grid on
box on
title('Width (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s2=nexttile(2);

momentsSpec.lrwidth(isnan(momentsSpec.lrwidth))=-99;
momentsSpec.lrwidth(isinf(momentsSpec.lrwidth))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lrwidth,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLrwidth);
s2.Colormap=col1;
colorbar
grid on
box on
title('Left edge to right edge width (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s4=nexttile(4);

momentsSpec.skew(isnan(momentsSpec.skew))=-99;
momentsSpec.skew(isinf(momentsSpec.skew))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s4.Colormap=col2;
colorbar
grid on
box on
title('Skewness (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s5=nexttile(5);

momentsSpec.lslope(isnan(momentsSpec.lslope))=-99;
momentsSpec.lslope(isinf(momentsSpec.lslope))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLslope);
s5.Colormap=col1;
colorbar
grid on
box on
title('Left slope (dB s m^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s6=nexttile(6);

momentsSpec.rslope(isnan(momentsSpec.rslope))=99;
momentsSpec.rslope(isinf(momentsSpec.rslope))=99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.rslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsRslope);
s6.Colormap=col1r;
colorbar
grid on
box on
title('Right slope (dB s m^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s7=nexttile(7);

momentsSpec.kurt(isnan(momentsSpec.kurt))=-99;
momentsSpec.kurt(isinf(momentsSpec.kurt))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s7.Colormap=col2;
colorbar
grid on
box on
title('Kurtosis (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s8=nexttile(8);

momentsSpec.level(isnan(momentsSpec.level))=-99;
momentsSpec.level(isinf(momentsSpec.level))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.level,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s8.Colormap=col2;
colorbar
grid on
box on
title('Left edge velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s9=nexttile(9);

momentsSpec.revel(isnan(momentsSpec.revel))=-99;
momentsSpec.revel(isinf(momentsSpec.revel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.revel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s9.Colormap=col2;
colorbar
grid on
box on
title('Right edge velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s10=nexttile(10);

momentsSpec.velRaw(isnan(momentsSpec.velRaw))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s10.Colormap=col2;
colorbar
grid on
box on
title('Velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s11=nexttile(11);

momentsSpec.lpvel(isnan(momentsSpec.lpvel))=-99;
momentsSpec.lpvel(isinf(momentsSpec.lpvel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.lpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s11.Colormap=col2;
colorbar
grid on
box on
title('Left peak velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s12=nexttile(12);

momentsSpec.rpvel(isnan(momentsSpec.rpvel))=-99;
momentsSpec.rpvel(isinf(momentsSpec.rpvel))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.rpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s12.Colormap=col2;
colorbar
grid on
box on
title('Right peak velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_specMomentsParams_',datestr(momentsSpec.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momentsSpec.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end