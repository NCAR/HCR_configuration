% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.1';
freqData='combined'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
ylimits=[0 6];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);
figdir=[indir(1:end-5)];

startTime=datetime(2018,1,22,21,40,0);
endTime=datetime(2018,1,22,21,56,0);

%% Load data

disp('Loading data ...');

data=[];

data.HSRL_Aerosol_Backscatter_Coefficient=[];

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
fig1=figure('DefaultAxesFontSize',13,'position',[100,100,1300,400]);

s1=subplot(1,1,1);
colormap(jet);

hold on;
sub1=surf(data.time,data.asl/1000,log10(data.HSRL_Aerosol_Backscatter_Coefficient),'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Aerosol backscatter coefficient (m^{-1} sr^{-1})')
grid on
box on
cb1=colorbar;
caxis([-8 -3])

s1.Position=[0.04,0.14,0.9,0.77];
cb1.Position=[0.95,0.14,0.0164,0.77];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'.png'],'-dpng','-r0');
