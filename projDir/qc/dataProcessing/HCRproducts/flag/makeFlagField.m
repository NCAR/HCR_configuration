% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2021,5,29,20,30,0);
endTime=datetime(2021,5,29,21,30,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset
quality='qc0'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v0.1';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc0/cfradial/v0.1/flagPlots/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

directories.dataDir=HCRdir(project,quality,qcVersion,freqData);

% if strcmp(whichModel,'era5')
%     directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',freqData,'/'];
% elseif strcmp(whichModel,'ecmwf')
%     directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/forecastInterp/',project,'/',freqData,'/'];
% end

%% Load data

data.DBZ=[];
%data.VEL=[];
%data.VEL_RAW=[];
%data.VEL_CORR=[];
data.WIDTH=[];
%data.WIDTH_CORR=[];
data.DBMVC=[];
%data.SNR=[];
%data.NCP=[];
%data.LDR=[];
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

%% Mask data

[maskData antStat]=echoMask(data);

refl=data.DBZ;
refl(maskData>1)=nan;

maskPlot=maskData;
maskPlot(maskData==0)=nan;

data.FLAG=maskData;

%% Plot
ylimits=[-0.5 15];

ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
    'Bang (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
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

figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200],'renderer','Zbuffer');

ax1=subplot(3,1,1);
fig1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
fig1=colMapDBZ(fig1);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity (all data)')

ax3=subplot(3,1,3);
fig3=surf(data.time,data.asl./1000,refl,'edgecolor','none');
view(2);
fig3=colMapDBZ(fig3);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity (cloud only, i.e. DBZ(FLAG>1)=NAN)')

ax2=subplot(3,1,2);
fig2=surf(data.time,data.asl./1000,maskPlot,'edgecolor','none');
view(2);
caxis([1 12]);
colormap(ax2,colMask);
hcb=colorbar;
set(hcb,'ytick',[1.5:0.91:12.5]);
set(hcb,'YTickLabel',ytickLabels);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Flag field')

linkaxes([ax1,ax2,ax3],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_echoMask'],...
   '-dpng','-r0');

%% Plot antenna status
figure('DefaultAxesFontSize',11,'position',[1,100,1800,300],'renderer','painters');

plot(data.time,antStat,'linewidth',2)
xlim([data.time(1),data.time(end)]);
title('Antenna status')
yticks(0:4)
yticklabels({'Down (0)','Up (1)','Pointing (2)','Scanning (3)','Transition (4)'})

set(gcf,'PaperPositionMode','auto')
print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_antStat'],...
   '-dpng','-r0');