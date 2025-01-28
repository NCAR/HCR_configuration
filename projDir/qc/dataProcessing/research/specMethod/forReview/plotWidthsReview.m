function plotWidthsReview(momentsSpecParams,momentsTimeW,momentsTimeSmoothW,momentsTimeWC,momentsTimeSmoothWC,figdir,project,showPlot)

aslGood=momentsSpecParams.asl(~isnan(momentsSpecParams.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,2];
climsDiff=[-1.5,1.5];
colWidth=cat(1,[0,0,0],jet);
colDiff=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsTimeW(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=colWidth;
colorbar
grid on
box on
title('Time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s2=nexttile(2);

momentsTimeWC(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeWC,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s2.Colormap=colWidth;
colorbar
grid on
box on
title('Time domain corr (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s3=nexttile(3);

momentsTimeSmoothW(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeSmoothW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s3.Colormap=colWidth;
colorbar
grid on
box on
title('Time domain smooth (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s4=nexttile(4);

momentsTimeSmoothWC(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeSmoothWC,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=colWidth;
colorbar
grid on
box on
title('Time domain smooth corr (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s6=nexttile(6);

momentsSpecParams.width(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsSpecParams.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s6.Colormap=colWidth;
colorbar
grid on
box on
title('Spectra corr (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'widthTest/',project,'/',project,'_widths_',datestr(momentsSpecParams.time(1),'yyyymmdd_HHMMSS'), ...
    '_to_',datestr(momentsSpecParams.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

%% Figure
f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsTimeW(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeW-momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain - time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s2=nexttile(2);

momentsTimeWC(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeWC-momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain corr - time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s3=nexttile(3);

momentsTimeSmoothW(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeSmoothW-momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain smooth - time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s4=nexttile(4);

momentsTimeSmoothWC(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsTimeSmoothWC-momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain smooth corr - time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

s6=nexttile(6);

momentsSpecParams.width(isnan(momentsSpecParams.velRaw))=-99;

hold on
surf(momentsSpecParams.time,momentsSpecParams.asl./1000,momentsSpecParams.width-momentsTimeW,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('Spectra corr - time domain (m s^{-1})')
ylim(ylims);
xlim([momentsSpecParams.time(1),momentsSpecParams.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'widthTest/',project,'/',project,'_widthDiffs_',datestr(momentsSpecParams.time(1),'yyyymmdd_HHMMSS'), ...
    '_to_',datestr(momentsSpecParams.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end