function plotWidths(moments,momentsSpBasic,momentsSpBasicRMnoise,momentsSpSmooth,momentsSpSmoothCorr,cf,plotTimeAll,figdir,project,showPlot)

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,2];
colDiff=cat(1,[0,0,0],jet);

%% Figure
f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

moments.width(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

% plotRangeInds=[20:20:700];
% for kk=1:length(plotTimeAll)
%     times=repmat(plotTimeAll(kk),length(plotRangeInds),1);
%     alts=momentsSpBasic.asl(plotRangeInds,momentsSpBasic.time==plotTimeAll(kk));
%     scatter(times,alts./1000,36,'k','+');
% end
% 
% s1.SortMethod='childorder';

s2=nexttile(2);

moments.widthCorr(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.widthCorr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpBasic.width(isnan(momentsSpBasic.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

momentsSpBasicRMnoise.width(isnan(momentsSpBasicRMnoise.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasicRMnoise.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s5=nexttile(5);

momentsSpSmooth.width(isnan(momentsSpSmooth.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmooth.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s5.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s6=nexttile(6);

momentsSpSmoothCorr.width(isnan(momentsSpSmoothCorr.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmoothCorr.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_width_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end