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
data.VEL_unfold_short=[];
data.WIDTH_short=[];
data.LDRV_short=[];
data.FLAG_short=[];
data.DBZ_long=[];
data.VEL_unfold_long=[];
data.WIDTH_long=[];
data.LDRV_long=[];
data.FLAG_long=[];
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

data.DBZ_short(data.FLAG_short~=1)=nan;
data.VEL_unfold_short(data.FLAG_short~=1)=nan;
data.WIDTH_short(data.FLAG_short~=1)=nan;
data.LDRV_short(data.FLAG_short~=1)=nan;

data.DBZ_long(data.FLAG_long~=1)=nan;
data.VEL_unfold_long(data.FLAG_long~=1)=nan;
data.WIDTH_long(data.FLAG_long~=1)=nan;
data.LDRV_long(data.FLAG_long~=1)=nan;

%% Plot DBZ and LDR
pix=50;
ylimDU=[0.2,10.5];

close all

f1 = figure('Position',[200 500 900 1000],'DefaultAxesFontSize',12);
t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight','TileIndexing', 'columnmajor');

%%% DBZ
s1=nexttile(1);
climVar=[-25,25];

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
colM=jet;
s1.Colormap=colM;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(a) DBZ_short','Interpreter','none');

s2=nexttile(2);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s2.Colormap=colM;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(c) DBZ_long','Interpreter','none');

s3=nexttile(3);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s3.Colormap=colM;
hcb=colorbar('Location','southoutside');
hcb.Title.String="dBZ";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

ylabel('Range (km)')
grid off
box on

title('(e) DBZ','Interpreter','none');

%%% LDR
s4=nexttile(4);
climVar=[-30,-15];
colL=jet;

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDRV_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s4.Colormap=colL;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

%ylabel('Range (km)')
grid off
box on

title('(b) LDR_short','Interpreter','none');

s5=nexttile(5);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDRV_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s5.Colormap=colL;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(d) LDR_long','Interpreter','none');

s6=nexttile(6);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.LDR(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s6.Colormap=colL;
hcb=colorbar('Location','southoutside');
hcb.Title.String="dB";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

%xticklabels('');
yticklabels('');

grid off
box on

title('(f) LDR','Interpreter','none');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'mergedExample1.png'],'-dpng','-r0')

%% Plot VEL and WIDTH
pix=50;
ylimDU=[0.2,10.5];

close all

f1 = figure('Position',[200 500 900 1000],'DefaultAxesFontSize',12);
t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight','TileIndexing', 'columnmajor');

%%% VEL
s1=nexttile(1);
climVar=[-10,10];

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL_unfold_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s1.Colormap=velCols;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');

ylabel('Range (km)')
grid off
box on

title('(a) VEL_short','Interpreter','none');

s2=nexttile(2);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL_unfold_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s2.Colormap=velCols;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);
ylabel('Range (km)')

xticklabels('');

grid off
box on

title('(c) VEL_long','Interpreter','none');

s3=nexttile(3);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.VEL(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s3.Colormap=velCols;
hcb=colorbar('Location','southoutside');
hcb.Title.String="m s^{-1}";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);
ylabel('Range (km)')

grid off
box on

title('(e) VEL','Interpreter','none');

%%% WIDTH
s4=nexttile(4);
climVar=[0,2];
colW=pink;

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH_short(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s4.Colormap=colW;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(b) WIDTH_short','Interpreter','none');

s5=nexttile(5);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s5.Colormap=colW;

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

xticklabels('');
yticklabels('');

grid off
box on

title('(d) WIDTH_long','Interpreter','none');

s6=nexttile(6);

surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.WIDTH(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim(climVar);
s6.Colormap=colW;
hcb=colorbar('Location','southoutside');
hcb.Title.String="m s^{-1}";

xlim([data.time(1),data.time(end)]);
ylim([ylimDU]);

yticklabels('');

grid off
box on

title('(f) WIDTH','Interpreter','none');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'mergedExample2.png'],'-dpng','-r0')