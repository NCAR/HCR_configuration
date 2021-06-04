% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2018,1,23,0,0,0);
endTime=datetime(2018,1,23,0,20,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.1';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

directories.dataDir=HCRdir(project,quality,qcVersion,freqData);

if strcmp(whichModel,'era5')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',freqData,'/'];
elseif strcmp(whichModel,'ecmwf')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/forecastInterp/',project,'/',freqData,'/'];
end

%% Load data

data.DBZ=[];
data.FLAG=[];
data.DBZ=[];
data.WIDTH=[];
data.DBMVC=[];
data.LDR=[];
data.TOPO=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Create FLAG
[maskData antStat]=echoMask(data);

refl=data.DBZ;
refl(maskData>1)=nan;

maskPlot=maskData;
maskPlot(maskData==0)=nan;

data.FLAG=maskData;
data.DBZ_MASKED=refl;

%% Plot
ylimits=[-0.5 6];

ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
    'Trans. Pulse (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
    'NS cal (10)';'Ant. trans. (11)';'Missing (12)'};

colMask=[0.4,0.8,1;
    0,0,0;
    0.5,0,0.5;
    0,1,0;
    0.2,0.6,0.2;
    1,0,0;
    0,0,0.6;
    0.7065,0.4415,0.2812;
    0.5,0.5,0.5;    
    0.9290,0.8940,0.1250;
    1,0,1;    
    1,0.6,0];

close all

wi=10;
hi=10;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1=subplot(3,1,1);
hold on;
outerpos1 = ax1.Position;
ax1.Position = [outerpos1(1)-0.07 outerpos1(2)+0.02 outerpos1(3)+0.04 outerpos1(4)+0.02];

fig1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
fig1=colMapDBZ(fig1);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(a)                      DBZ')

ax3=subplot(3,1,3);
hold on;
outerpos1 = ax3.Position;
ax3.Position = [outerpos1(1)-0.07 outerpos1(2)-0.04 outerpos1(3)+0.04 outerpos1(4)+0.02];
fig3=surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
view(2);
fig3=colMapDBZ(fig3);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(c)        DBZ_MASKED','interpreter','none')

ax2=subplot(3,1,2);
hold on;
outerpos1 = ax2.Position;
ax2.Position = [outerpos1(1)-0.07 outerpos1(2)-0.01 outerpos1(3)+0.04 outerpos1(4)+0.02];
fig2=surf(data.time,data.asl./1000,data.FLAG,'edgecolor','none');
view(2);
caxis([1 12]);
colormap(ax2,colMask);
hcb=colorbar;
set(hcb,'ytick',[1.5:0.91:12.5]);
set(hcb,'YTickLabel',ytickLabels);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(b)                    FLAG')

print([figdir,'maskExample.png'],'-dpng','-r0');