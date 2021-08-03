% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
if strcmp(project,'otrec')
    ylimits=[-0 15];
elseif strcmp(project,'socrates')
    ylimits=[-0.2 6];
elseif strcmp(project,'cset')
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/HCR/other/'];
%figdir=['/home/romatsch/plots/HCR/meltingLayer/selected/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,qcVersion,freqData);

startTime=datetime(2019,10,1,16,2,0);
endTime=datetime(2019,10,1,16,15,30);

%% Load data

disp('Loading data ...');

data=[];

data.DBZ=[];
data.FLAG=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

data.DBZ(data.FLAG>1)=nan;

%% Find melting layer
close all
fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);

colormap(jet);

hold on;
sub1=surf(data.time,data.asl/1000,data.DBZ,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity (dBZ)')
grid on
box on
colorbar
caxis([-50 20])

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
