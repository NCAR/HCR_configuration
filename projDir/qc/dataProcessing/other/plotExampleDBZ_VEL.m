% Analyze HCR clouds

clear all;
close all;

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
ylimits=[1 13.5];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/HCR/other/'];
%figdir=['/home/romatsch/plots/HCR/meltingLayer/selected/',project,'/'];


indir=HCRdir(project,quality,qcVersion,freqData);

startTime=datetime(2021,6,21,1,6,0);
endTime=datetime(2021,6,21,1,21,0);

%% Load data

disp('Loading data ...');

data=[];

data.DBZ_MASKED=[];
data.VEL_MASKED=[];

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

%% Plot

close all
fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,800]);

s1=subplot(2,1,1);
colormap(jet);

hold on;
sub1=surf(data.time,data.asl/1000,data.DBZ_MASKED,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity (dBZ)')
grid on
box on
colorbar
caxis([-40 30])

s2=subplot(2,1,2);
s2.Colormap=flipud(jet);

hold on;
sub1=surf(data.time,data.asl/1000,data.VEL_MASKED,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Radial velocity (m s^{-1})')
grid on
box on
colorbar
caxis([-12 12])

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
