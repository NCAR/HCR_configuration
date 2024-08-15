function plotKurtosis(moments,momentsSpBasic,momentsSpBasicRMnoise,momentsSpSmooth,momentsSpSmoothCorr,cf,plotTimeAll,figdir,project,showPlot)

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.kurt))./1000;
ylims=[0,max(aslGood)+0.5];

climsKurt=[-6,6];
colDiff=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1600 830],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsSpBasic.kurt(isnan(momentsSpBasic.kurt))=-99;
momentsSpBasic.kurt(isinf(momentsSpBasic.kurt))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s2=nexttile(2);

momentsSpBasicRMnoise.kurt(isnan(momentsSpBasicRMnoise.kurt))=-99;
momentsSpBasicRMnoise.kurt(isinf(momentsSpBasicRMnoise.kurt))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasicRMnoise.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpSmooth.kurt(isnan(momentsSpSmooth.kurt))=-99;
momentsSpSmooth.kurt(isinf(momentsSpSmooth.kurt))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmooth.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

momentsSpSmoothCorr.kurt(isnan(momentsSpSmoothCorr.kurt))=-99;
momentsSpSmoothCorr.kurt(isinf(momentsSpSmoothCorr.kurt))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmoothCorr.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_kurt_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end