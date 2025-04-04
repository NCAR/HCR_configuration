clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';

indir=HCRdir(project,quality,qcVersion,freqData);

%% Run processing


startTime=datetime(2024,6,14,20,27,0);
endTime=datetime(2024,6,14,21,32,0);

data=[];

data.SNRVC_short=[];
data.DBZ=[];

%% Load data
disp('Loading data ...');

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Load other
load([figdir,'mergeScatterCensor.mat']);

edges.DBZ=[-100,-40:1:20,50];
edges.VEL=[-30,-12:0.25:12,30];
edges.WIDTH=[0:0.1:4,10];
edges.LDRV=[-50,-30:1:10,30];

%% Plot
close all
f1 = figure('Position',[200 500 500 600],'DefaultAxesFontSize',8);
t = tiledlayout(3,2,'TileSpacing','tight','Padding','compact');
col=cat(1,[1,1,1],jet);

thisName='DBZ';
[N,~,~]=histcounts2(bySNRall.bin1.(thisName)(:,1), ...
    bySNRall.bin1.(thisName)(:,2),edges.(thisName),edges.(thisName));

plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

s1=nexttile(1);
hold on
s1.Colormap=col;
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
s2.Colormap=col;
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
xlabel('VEL_long','Interpreter','none')
ylabel('VEL_short','Interpreter','none')
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
s3.Colormap=col;
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
s4.Colormap=col;
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

%%% SNR
pix=50;
ylimDU=[0.2,10.5];

s5=nexttile(6);

hold on
surf(1:pix:length(data.time),data.range(:,1)./1000,data.SNRVC_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim([1000,1100]);
colM=colormap('cool');
colormap(colM);

SNRcont=data.SNRVC_short;
SNRcont(isnan(data.DBZ))=nan;
%ylabel('Range (km)')
yticklabels('');

%xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

yyaxis right
contourf(1:length(data.time),data.range(:,1)./1000,SNRcont,[-5.5,10,25],'LineColor','none')
s5.SortMethod='childorder';
xlim([1,length(data.time)]);
yticks('')
grid off
box on
colCont=cool(3);
colCont=[[0,0,0];colCont];
s5.Colormap=colCont;
set(gca,'YColor','k');
ylim([ylimDU]);

l1=plot(nan,'-','Color',[0,1,1],'LineWidth',2);
l2=plot(nan,'-','Color',[0.5,0.5,1],'LineWidth',2);
l3=plot(nan,'-','Color',[1,0,1],'LineWidth',2);

l=legend([l1,l2,l3],{'DBZ/VEL','WIDTH','LDR'});
l.IconColumnWidth = 10;
l.Location='southoutside';
l.Orientation='horizontal';

title('(f) SNR thresholds','Interpreter','none');
yticklabels('');
xticklabels('');

s6=nexttile(5);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.SNRVC_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim([-10,55]);
colM=jet;
s6.Colormap=colM;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

ylabel('Range (km)')
grid off
box on

title('(e) SNR_short','Interpreter','none');
hcb=colorbar;
hcb.Title.String="dB";

s1.Colormap=col;
s2.Colormap=col;
s3.Colormap=col;
s4.Colormap=col;

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'mergedScatter.png'],'-dpng','-r0')