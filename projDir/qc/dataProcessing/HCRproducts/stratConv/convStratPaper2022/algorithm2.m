% Plot HCR convStrat from mat file in hourly plots

clear all;
close all;

ylimUpper=14.99;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

colmapBasic=[0,0,1;
    0.32,0.78,0.59;
    0.7,0,0];

load('/scr/sci/romatsch/other/convStratPaperHCR/algPlotData.mat');

%% Plot

close all

wi=9;
hi=4.5;

fig1=figure('DefaultAxesFontSize',13,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(3,1,1);

surf(data.time,data.asl./1000,classBasic,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
s1.Colormap=colmapBasic;
caxis([0.5 3.5]);
cb1=colorbar;
cb1.Ticks=1:3;
cb1.TickLabels={'Stratiform','Mixed','Convective'};
set(gca,'XTickLabel',[]);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(a) Basic echo type','FontSize',11,'FontWeight','bold');

s2=subplot(3,1,2);
hold on
surf(data.time,data.asl./1000,classSubPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
s2.Colormap=colmapSC;
caxis([0.5 9.5]);
cb2=colorbar;
cb2.Ticks=1:9;
cb2.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(b) Advanced echo type','FontSize',11,'FontWeight','bold');

s3=subplot(3,1,3);

hold on
scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
%set(gca,'clim',[0,1]);
set(gca,'YTickLabel',[]);
s3.Colormap=colmapSC;
xlim([data.time(1),data.time(end)]);
grid on
box on

s1.Position=[0.06 0.58 0.785 0.4];
s2.Position=[0.06 0.15 0.785 0.4];
s3.Position=[0.06 0.09 0.785 0.03];

cb1.Position=[0.855,0.58,0.02,0.4];
cb2.Position=[0.855,0.15,0.02,0.4];


set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'algorithm2.tif'],'-dtiffn','-r0')
