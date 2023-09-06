% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
ylimits=[0 14];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);
figdir=[indir(1:end-5)];

startTime=datetime(2019,10,1,14,35,0);
endTime=datetime(2019,10,1,14,50,0);

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
disp('Plotting ...');

close all
fig1=figure('DefaultAxesFontSize',13,'position',[100,100,1300,800]);

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
cb1=colorbar;
caxis([-20 15])

s2=subplot(2,1,2);

hold on;
sub2=surf(data.time,data.asl/1000,data.VEL_MASKED,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Radial velocity (m s^{-1})')
grid on
box on
s2.Colormap=velCols;
caxis([-5 5]);
cb2=colorbar;
cb2.Ticks=-5:5;

s1.Position=[0.04,0.56,0.9,0.4];
s2.Position=[0.04,0.07,0.9,0.4];

cb1.Position=[0.95,0.56,0.0164,0.4];
cb2.Position=[0.95,0.07,0.0164,0.4];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
