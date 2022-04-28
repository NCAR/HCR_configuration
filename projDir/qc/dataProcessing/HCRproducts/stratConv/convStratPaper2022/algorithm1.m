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
hi=10;

fig1=figure('DefaultAxesFontSize',13,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(5,1,1);

colormap jet

hold on
surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
set(gca,'XTickLabel',[]);
cb1=colorbar;
plot(data.time,data.altitude./1000,'-k','LineWidth',2);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

s2=subplot(5,1,2);

colormap jet

hold on
surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-5 5]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
set(gca,'XTickLabel',[]);
cb2=colorbar;
cb2.Ticks=-4:2:4;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(b) Velocity (m s^{-1})','FontSize',11,'FontWeight','bold');

s3=subplot(5,1,3);

colormap jet

hold on
surf(data.time,data.asl./1000,convDBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
set(gca,'XTickLabel',[]);
cb3=colorbar;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(c) Scaled TDBZ','FontSize',11,'FontWeight','bold');

s4=subplot(5,1,4);

hold on
surf(data.time,data.asl./1000,convVEL,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
caxis([0 1]);
set(gca,'XTickLabel',[]);
cb4=colorbar;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(d) Scaled TVEL','FontSize',11,'FontWeight','bold');

s5=subplot(5,1,5);

hold on
surf(data.time,data.asl./1000,convBoth,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
caxis([0 1]);
cb5=colorbar;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(e) Convectivity','FontSize',11,'FontWeight','bold');

s1.Position=[0.06 0.81 0.87 0.18];
s2.Position=[0.06 0.62 0.87 0.18];
s3.Position=[0.06 0.43 0.87 0.18];
s4.Position=[0.06 0.24 0.87 0.18];
s5.Position=[0.06 0.05 0.87 0.18];

cb1.Position=[0.94,0.815,0.02,0.17];
cb2.Position=[0.94,0.625,0.02,0.17];
cb3.Position=[0.94,0.435,0.02,0.17];
cb4.Position=[0.94,0.245,0.02,0.17];
cb5.Position=[0.94,0.055,0.02,0.17];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'algorithm1.tif'],'-dtiffn','-r0')
