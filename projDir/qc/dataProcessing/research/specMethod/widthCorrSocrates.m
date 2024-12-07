clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'widthSocrates.mat']);

%% Figure

close all

aslGood=momentsSpecBasic.asl(~isnan(momentsSpecBasic.velRaw))./1000;
ylims=[0,2.2];
xlims=[datetime(2018,2,4,2,47,15),datetime(2018,2,4,2,49,5)];

climsWidth=[0,1.2];
climsDbz=[-15,25];
colDiff=jet;

% Figure
f1 = figure('Position',[200 500 1000 600],'DefaultAxesFontSize',12);

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsTime.widthCorr(isnan(dataCF.VEL_MASKED))=nan;

hold on
surf(momentsTime.time,momentsTime.asl./1000,momentsTime.widthCorr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('(a) WIDTH pulse-pair (m s^{-1})')
ylim(ylims);
xlim(xlims);

s2=nexttile(2);

momentsSpecSmoothCorr.width(isnan(momentsSpecSmoothCorr.width))=nan;

hold on
surf(momentsSpecBasic.time,momentsSpecBasic.asl./1000,momentsSpecSmoothCorr.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('(b) WIDTH spectra (m s^{-1})')
ylim(ylims);
xlim(xlims);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'widthCorr_socrates.png'],'-dpng','-r0');
