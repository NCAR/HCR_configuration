clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

load([figdir,'cfad.mat']);

% Long
f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12);
t = tiledlayout(2,5,'TileSpacing','tight','Padding','compact');
colormap('jet');

surf(dbzEdges(1:end-1),range(:,1)./1000,countsIOPallL,'edgecolor','none');
view(2);

clim([0,20000])

xlim([-50,20])

colorbar

xlabel('Reflectivity (dBZ)')
ylabel('Range (km)')

box on
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'cfad'],'-dpng','-r0')
