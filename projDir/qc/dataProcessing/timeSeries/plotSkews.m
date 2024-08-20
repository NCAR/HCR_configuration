function plotSkews(moments,momentsSpBasic,momentsSpBasicRMnoise,momentsSpSmooth,momentsSpSmoothCorr,cf,plotTimeAll,figdir,project,showPlot)
% momentsSpBasic.skew(:,momentsSpBasic.elevation>0,:)=-momentsSpBasic.skew(:,momentsSpBasic.elevation>0,:);
% momentsSpBasicRMnoise.skew(:,momentsSpBasic.elevation>0,:)=-momentsSpBasicRMnoise.skew(:,momentsSpBasic.elevation>0,:);
% momentsSpSmooth.skew(:,momentsSpBasic.elevation>0,:)=-momentsSpSmooth.skew(:,momentsSpBasic.elevation>0,:);
% momentsSpSmoothCorr.skew(:,momentsSpBasic.elevation>0,:)=-momentsSpSmoothCorr.skew(:,momentsSpBasic.elevation>0,:);

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.skew))./1000;
ylims=[0,max(aslGood)+0.5];

climsSkew=[-3,3];
colDiff=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1600 830],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsSpBasic.skew(isnan(momentsSpBasic.skew))=-99;
momentsSpBasic.skew(isinf(momentsSpBasic.skew))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s2=nexttile(2);

momentsSpBasicRMnoise.skew(isnan(momentsSpBasicRMnoise.skew))=-99;
momentsSpBasicRMnoise.skew(isinf(momentsSpBasicRMnoise.skew))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasicRMnoise.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpSmooth.skew(isnan(momentsSpSmooth.skew))=-99;
momentsSpSmooth.skew(isinf(momentsSpSmooth.skew))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmooth.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

momentsSpSmoothCorr.skew(isnan(momentsSpSmoothCorr.skew))=-99;
momentsSpSmoothCorr.skew(isinf(momentsSpSmoothCorr.skew))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmoothCorr.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_skew_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end