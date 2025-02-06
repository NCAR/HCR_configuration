clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'highMomentExample_spicule.mat']);

momentsSpecSmoothCorr.velRaw(:,momentsSpecSmoothCorr.elevation>0)=-momentsSpecSmoothCorr.velRaw(:,momentsSpecSmoothCorr.elevation>0);
momentsSpecSmoothCorr.skew(:,momentsSpecSmoothCorr.elevation>0)=-momentsSpecSmoothCorr.skew(:,momentsSpecSmoothCorr.elevation>0);

aslGood=momentsSpecSmoothCorr.asl(~isnan(momentsSpecSmoothCorr.velRaw))./1000;
xlims=[datetime(2021,5,29,16,19,45),datetime(2021,5,29,16,23,53)];
ylims=[3.3,7.8];

climsLrwidth=[0,12];
climsLslope=[0,18];
climsRslope=[-18,0];
climsVel=[-10,10];

lslopeOrig=momentsSpecSmoothCorr.lslope;
rslopeOrig=momentsSpecSmoothCorr.rslope;
momentsSpecSmoothCorr.lslope(:,momentsSpecSmoothCorr.elevation>0)=-rslopeOrig(:,momentsSpecSmoothCorr.elevation>0);
momentsSpecSmoothCorr.rslope(:,momentsSpecSmoothCorr.elevation>0)=-lslopeOrig(:,momentsSpecSmoothCorr.elevation>0);

levelOrig=momentsSpecSmoothCorr.level;
revelOrig=momentsSpecSmoothCorr.revel;
momentsSpecSmoothCorr.level(:,momentsSpecSmoothCorr.elevation>0)=-revelOrig(:,momentsSpecSmoothCorr.elevation>0);
momentsSpecSmoothCorr.revel(:,momentsSpecSmoothCorr.elevation>0)=-levelOrig(:,momentsSpecSmoothCorr.elevation>0);

lpvelOrig=momentsSpecSmoothCorr.lpvel;
rpvelOrig=momentsSpecSmoothCorr.rpvel;
momentsSpecSmoothCorr.lpvel(:,momentsSpecSmoothCorr.elevation>0)=-rpvelOrig(:,momentsSpecSmoothCorr.elevation>0);
momentsSpecSmoothCorr.rpvel(:,momentsSpecSmoothCorr.elevation>0)=-lpvelOrig(:,momentsSpecSmoothCorr.elevation>0);

momentsSpecSmoothCorr.lpvel(isnan(momentsSpecSmoothCorr.rpvel))=nan;
momentsSpecSmoothCorr.rpvel(isnan(momentsSpecSmoothCorr.lpvel))=nan;

col1=jet;
col1r=flipud(col1);
col2=velCols;

%% Figure
f1 = figure('Position',[200 500 1400 1200],'DefaultAxesFontSize',12,'visible','on');

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

momentsSpecSmoothCorr.lrwidth(isnan(momentsSpecSmoothCorr.lrwidth))=nan;
momentsSpecSmoothCorr.lrwidth(isinf(momentsSpecSmoothCorr.lrwidth))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.lrwidth,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLrwidth);
s1.Colormap=col1;
colorbar
grid on
box on
title('(a) Edge-to-edge width (m s^{-1})')
ylim(ylims);
xlim(xlims);

s3=nexttile(3);

momentsSpecSmoothCorr.lslope(isnan(momentsSpecSmoothCorr.lslope))=nan;
momentsSpecSmoothCorr.lslope(isinf(momentsSpecSmoothCorr.lslope))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.lslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsLslope);
s3.Colormap=col1;
colorbar
grid on
box on
title('(b) Left slope (dB s m^{-1})')
ylim(ylims);
xlim(xlims);

s4=nexttile(4);

momentsSpecSmoothCorr.rslope(isnan(momentsSpecSmoothCorr.rslope))=nan;
momentsSpecSmoothCorr.rslope(isinf(momentsSpecSmoothCorr.rslope))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.rslope,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsRslope);
s4.Colormap=col1r;
colorbar
grid on
box on
title('(c) Right slope (dB s m^{-1})')
ylim(ylims);
xlim(xlims);

s5=nexttile(5);

momentsSpecSmoothCorr.level(isnan(momentsSpecSmoothCorr.level))=nan;
momentsSpecSmoothCorr.level(isinf(momentsSpecSmoothCorr.level))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.level,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s5.Colormap=col2;
colorbar
grid on
box on
title('(d) Left-edge velocity (m s^{-1})')
ylim(ylims);
xlim(xlims);

s6=nexttile(6);

momentsSpecSmoothCorr.revel(isnan(momentsSpecSmoothCorr.revel))=nan;
momentsSpecSmoothCorr.revel(isinf(momentsSpecSmoothCorr.revel))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.revel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s6.Colormap=col2;
colorbar
grid on
box on
title('(e) Right-edge velocity (m s^{-1})')
ylim(ylims);
xlim(xlims);

s7=nexttile(7);

momentsSpecSmoothCorr.lpvel(isnan(momentsSpecSmoothCorr.lpvel))=nan;
momentsSpecSmoothCorr.lpvel(isinf(momentsSpecSmoothCorr.lpvel))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.lpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s7.Colormap=col2;
colorbar
grid on
box on
title('(f) Left-peak velocity (m s^{-1})')
ylim(ylims);
xlim(xlims);

s8=nexttile(8);

momentsSpecSmoothCorr.rpvel(isnan(momentsSpecSmoothCorr.rpvel))=nan;
momentsSpecSmoothCorr.rpvel(isinf(momentsSpecSmoothCorr.rpvel))=nan;

hold on
surf(momentsSpecSmoothCorr.time,momentsSpecSmoothCorr.asl./1000,momentsSpecSmoothCorr.rpvel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s8.Colormap=col2;
colorbar
grid on
box on
title('(g) Right-peak velocity (m s^{-1})')
ylim(ylims);
xlim(xlims);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'specParamsSPICULE.png'],'-dpng','-r0');