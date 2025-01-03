clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

load([figdir,'mergeScatter.mat']);

edges.DBZ=[-100,-40:1:20,50];
edges.VEL=[-30,-12:0.25:12,30];
edges.WIDTH=[0:0.1:4,10];
edges.LDRV=[-50,-30:1:10,30];

%% Plot
close all
f1 = figure('Position',[200 500 800 600],'DefaultAxesFontSize',12);
t = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
col=cat(1,[1,1,1],jet);
colormap(col);

thisName='DBZ';
[N,~,~]=histcounts2(bySNRall.bin1.(thisName)(:,1), ...
    bySNRall.bin1.(thisName)(:,2),edges.(thisName),edges.(thisName));

plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

s1=nexttile(1);
hold on
pcolor(plotCoords,plotCoords,log(N));
shading('flat');
caxis([0,12])
floorAx=-39;
ceilAx=-6;
xlim([floorAx,ceilAx]);
ylim([floorAx,ceilAx]);
plot([floorAx,ceilAx],[floorAx,ceilAx],'-k','LineWidth',2);
box on
%grid on
xlabel('DBZ_long','Interpreter','none')
ylabel('DBZ_short','Interpreter','none')
s1.SortMethod='childorder';
set(gca,'layer','top')
title('(a) DBZ, -6 to -5 dB SNR');
hcb=colorbar;
hcb.Title.String="log";
xticks(-100:5:100);
yticks(-100:5:100);

thisName='VEL';
[N,~,~]=histcounts2(bySNRall.bin1.(thisName)(:,1), ...
    bySNRall.bin1.(thisName)(:,2),edges.(thisName),edges.(thisName));

plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

s2=nexttile(2);
hold on
pcolor(plotCoords,plotCoords,log(N));
shading('flat');
caxis([0,12])
floorAx=-5;
ceilAx=5;
xlim([floorAx,ceilAx]);
ylim([floorAx,ceilAx]);
plot([floorAx,ceilAx],[floorAx,ceilAx],'-k','LineWidth',2);
box on
%grid on
xlabel('DBZ_long','Interpreter','none')
ylabel('DBZ_short','Interpreter','none')
s2.SortMethod='childorder';
set(gca,'layer','top')
title('(b) VEL, -5 to -4 dB SNR');
hcb=colorbar;
hcb.Title.String="log";
xticks(-10:2:10);
yticks(-10:2:10);

thisName='WIDTH';
[N,~,~]=histcounts2(bySNRall.bin1.(thisName)(:,1), ...
    bySNRall.bin1.(thisName)(:,2),edges.(thisName),edges.(thisName));

plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

s3=nexttile(3);
hold on
pcolor(plotCoords,plotCoords,log(N));
shading('flat');
caxis([0,12])
floorAx=0.2;
ceilAx=3;
xlim([floorAx,ceilAx]);
ylim([floorAx,ceilAx]);
plot([floorAx,ceilAx],[floorAx,ceilAx],'-k','LineWidth',2);
box on
%grid on
xlabel('WIDTH_long','Interpreter','none')
ylabel('WIDTH_short','Interpreter','none')
s3.SortMethod='childorder';
set(gca,'layer','top')
title('(c) WIDTH, 10 to 11 dB SNR');
hcb=colorbar;
hcb.Title.String="log";
xticks(-10:0.5:10);
yticks(-10:0.5:10);

thisName='LDRV';
[N,~,~]=histcounts2(bySNRall.bin1.(thisName)(:,1), ...
    bySNRall.bin1.(thisName)(:,2),edges.(thisName),edges.(thisName));

plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

s4=nexttile(4);
hold on
pcolor(plotCoords,plotCoords,log(N));
shading('flat');
caxis([0,12])
floorAx=-28;
ceilAx=-5;
xlim([floorAx,ceilAx]);
ylim([floorAx,ceilAx]);
plot([floorAx,ceilAx],[floorAx,ceilAx],'-k','LineWidth',2);
box on
%grid on
xlabel('LDR_long','Interpreter','none')
ylabel('LDR_short','Interpreter','none')
s4.SortMethod='childorder';
set(gca,'layer','top')
title('(d) LDR, 25 to 26 dB SNR');
hcb=colorbar;
hcb.Title.String="log";
xticks(-100:5:100);
yticks(-100:5:100);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'mergedScatter.png'],'-dpng','-r0')