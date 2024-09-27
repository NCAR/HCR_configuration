clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'widthNoreaster.mat']);

%% Figure

close all

aslGood=momentsSpecBasic.asl(~isnan(momentsSpecBasic.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,1.7];
climsDbz=[-15,25];
colDiff=jet;

% Figure
f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12);

t = tiledlayout(4,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsTime.width(isnan(dataCF.VEL_MASKED))=nan;

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.width,'edgecolor','none');
view(2);
l=plot(dataCF.time,dataCF.altitude./1000,'-b','LineWidth',2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('(a) WIDTH_{raw} pulse-pair (m s^{-1})')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);
legend(l,'Aircraft altitude','Location','southeast');

s2=nexttile(2);

momentsTime.widthCorr(isnan(dataCF.VEL_MASKED))=nan;

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.widthCorr,'edgecolor','none');
view(2);
l=plot(dataCF.time,dataCF.altitude./1000,'-b','LineWidth',2);
ylabel('Altitude (km)');
clim(climsWidth);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('(b) WIDTH pulse-pair (m s^{-1})')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);
legend(l,'Aircraft altitude','Location','southeast');

s3=nexttile(3);

momentsSpecSmoothCorr.width(isnan(momentsSpecSmoothCorr.width))=nan;

hold on
surf(momentsSpecBasic.time,momentsSpecBasic.asl./1000,momentsSpecSmoothCorr.width,'edgecolor','none');
view(2);
l=plot(dataCF.time,dataCF.altitude./1000,'-b','LineWidth',2);
ylabel('Altitude (km)');
clim(climsWidth);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('(c) WIDTH spectra (m s^{-1})')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);
legend(l,'Aircraft altitude','Location','southeast');

s4=nexttile(4);

dataCF.DBZ(isnan(dataCF.VEL_MASKED))=nan;

hold on
surf(momentsTime.time,momentsTime.asl./1000,dataCF.DBZ,'edgecolor','none');
view(2);
l=plot(dataCF.time,dataCF.altitude./1000,'-b','LineWidth',2);
ylabel('Altitude (km)');
clim(climsDbz);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('(d) Reflectivity (dBZ)')
ylim(ylims);
xlim([momentsTime.time(1),momentsTime.time(end)]);
legend(l,'Aircraft altitude','Location','southeast');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'widthCorr_noreaster.png'],'-dpng','-r0');
