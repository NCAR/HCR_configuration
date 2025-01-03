% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';

indir=HCRdir(project,quality,qcVersion,freqData);

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

%% Run processing


startTime=datetime(2024,6,14,20,27,0);
endTime=datetime(2024,6,14,21,32,0);

data=[];

data.DBZ_short=[];
data.VEL_short=[];
data.WIDTH_short=[];
data.SNRVC_short=[];
data.LDRV_short=[];
data.DBZ_long=[];
data.VEL_long=[];
data.WIDTH_long=[];
data.LDRV_long=[];
data.DBZ=[];
data.VEL=[];
data.WIDTH=[];
data.LDR=[];

%% Load data
disp('Loading data ...');

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

%% Plot
pix=50;
ylimDU=[0.2,10.5];

close all

f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12);
t = tiledlayout(5,3,'TileSpacing','compact','Padding','tight');


%%% SNR
s2=nexttile(2);

hold on
surf(1:pix:length(data.time),data.range(:,1)./1000,data.SNRVC_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim([1000,1100]);
colM=colormap('cool');
colormap(colM);

SNRcont=data.SNRVC_short;
SNRcont(isnan(data.DBZ))=nan;

%xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

yyaxis right
contourf(1:length(data.time),data.range(:,1)./1000,SNRcont,[-5.5,10,25],'LineColor','none')
s2.SortMethod='childorder';
xlim([1,length(data.time)]);
%ylabel('Range (km)')
yticks('')
grid off
box on
colCont=cool(3);
colCont=[[0,0,0];colCont];
s2.Colormap=colCont;
set(gca,'YColor','k');
ylim([ylimDU]);

l1=plot(nan,'-','Color',[0,1,1],'LineWidth',2);
l2=plot(nan,'-','Color',[0.5,0.5,1],'LineWidth',2);
l3=plot(nan,'-','Color',[1,0,1],'LineWidth',2);

legend([l1,l2,l3],{'DBZ/VEL','WIDTH','LDR'});

title('(a) SNR thresholds','Interpreter','none');

s3=nexttile(3);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.SNRVC_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim([-10,55]);
colM=jet;
s3.Colormap=colM;
%hcb=colorbar;
%hcb.Title.String="dB";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

%ylabel('Range (km)')
grid off
box on

title('(b) SNR_short','Interpreter','none');
hcb=colorbar;
hcb.Title.String="dB";

%%% DBZ
s4=nexttile(4);
climVar=[-25,25];

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s4.Colormap=colM;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(c) DBZ_short','Interpreter','none');

s5=nexttile(5);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s5.Colormap=colM;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(d) DBZ_long','Interpreter','none');

s6=nexttile(6);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s6.Colormap=colM;
hcb=colorbar;
hcb.Title.String="dBZ";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(e) DBZ','Interpreter','none');


%%% VEL
s7=nexttile(7);
climVar=[-10,10];

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s7.Colormap=velCols;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(f) VEL_short','Interpreter','none');

s8=nexttile(8);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s8.Colormap=velCols;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(g) VEL_long','Interpreter','none');

s9=nexttile(9);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s9.Colormap=velCols;
hcb=colorbar;
hcb.Title.String="m s^{-1}";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(h) VEL','Interpreter','none');

%%% WIDTH
s10=nexttile(10);
climVar=[0,2];
colW=pink;

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s10.Colormap=colW;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(i) WIDTH_short','Interpreter','none');

s11=nexttile(11);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s11.Colormap=colW;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(j) WIDTH_long','Interpreter','none');

s12=nexttile(12);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s12.Colormap=colW;
hcb=colorbar;
hcb.Title.String="m s^{-1}";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(k) WIDTH','Interpreter','none');

%%% LDR
s13=nexttile(13);
climVar=[-30,-15];
colL=jet;

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDRV_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s13.Colormap=colL;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

%xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(l) LDR_short','Interpreter','none');

s14=nexttile(14);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDRV_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s14.Colormap=colL;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

%xticklabels('');
yticklabels('');

grid off
box on

title('(m) LDR_long','Interpreter','none');

s15=nexttile(15);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDR(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s15.Colormap=colL;
hcb=colorbar;
hcb.Title.String="dB";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

%xticklabels('');
yticklabels('');

grid off
box on

title('(n) LDR','Interpreter','none');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'mergedExample.png'],'-dpng','-r0')