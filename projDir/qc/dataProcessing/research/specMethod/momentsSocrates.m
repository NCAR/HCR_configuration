clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'highMomentExample_socrates.mat']);

momentsSpecSmoothCorr.velRaw(:,momentsSpecSmoothCorr.elevation>0)=-momentsSpecSmoothCorr.velRaw(:,momentsSpecSmoothCorr.elevation>0);
momentsSpecSmoothCorr.skew(:,momentsSpecSmoothCorr.elevation>0)=-momentsSpecSmoothCorr.skew(:,momentsSpecSmoothCorr.elevation>0);

aslGood=momentsSpecSmoothCorr.asl(~isnan(momentsSpecSmoothCorr.velRaw))./1000;
xlims=[datetime(2018,2,24,3,0,46),datetime(2018,2,24,3,2,0)];
ylims=[0,2.8];

climsDbz=[-15,25];
climsSnr=[-10,30];
climsVel=[-10,10];
climsWidth=[0,2.5];
climsSkew=[-3,3];
climsKurt=[-6,6];

col1=jet;
col2=velCols;

%% Figure
f1 = figure('Position',[200 500 1400 900],'DefaultAxesFontSize',12,'visible','on');

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

dataCF.DBZ(isnan(dataCF.VEL_MASKED))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,dataCF.DBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDbz);
s1.Colormap=col1;
colorbar
grid on
box on
title('(a) Reflectivity (dBZ)')
ylim(ylims);
xlim(xlims);

s2=nexttile(2);

dataCF.SNR(isnan(dataCF.VEL_MASKED))=nan;
dataCF.SNR(isnan(dataCF.SNR))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,dataCF.SNR,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSnr);
s2.Colormap=col1;
colorbar
grid on
box on
title('(b) Signal-to-noise ratio (dB)')
ylim(ylims);
xlim(xlims);

s3=nexttile(3);

momentsSpecSmoothCorr.velRaw(isnan(momentsSpecSmoothCorr.velRaw))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s3.Colormap=col2;
colorbar
grid on
box on
title('(c) Doppler velocity (m s^{-1})')
ylim(ylims);
xlim(xlims);

s4=nexttile(4);

momentsSpecSmoothCorr.width(isnan(momentsSpecSmoothCorr.width))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=col1;
colorbar
grid on
box on
title('(d) Spectrum width (m s^{-1})')
ylim(ylims);
xlim(xlims);

s5=nexttile(5);

momentsSpecSmoothCorr.skew(isnan(momentsSpecSmoothCorr.skew))=nan;
momentsSpecSmoothCorr.skew(isinf(momentsSpecSmoothCorr.skew))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s5.Colormap=col2;
colorbar
grid on
box on
title('(e) Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s6=nexttile(6);

momentsSpecSmoothCorr.kurt(isnan(momentsSpecSmoothCorr.kurt))=nan;
momentsSpecSmoothCorr.kurt(isinf(momentsSpecSmoothCorr.kurt))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s6.Colormap=col2;
colorbar
grid on
box on
title('(f) Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);


set(gcf,'PaperPositionMode','auto')
%print(f1,[figdir,'momentsSOCRATES.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'momentsSOCRATES.png'],'Resolution','300');