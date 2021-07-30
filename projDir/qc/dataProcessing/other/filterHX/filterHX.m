% Remove stripes in HX

clear all;
close all;

startTime=datetime(2021,6,11,18,45,0);
endTime=datetime(2021,6,11,19,5,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.0/flagPlots/'];

% if ~exist(figdir, 'dir')
%     mkdir(figdir)
% end

dataDir=HCRdir(project,quality,qcVersion,freqData);
%dataDir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/moments/10hz/';

%% Load data

data.DBMHX=[];
data.LDR=[];
data.FLAG=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

%% Find indices where correction needs to be applied
medianHXend=median(data.DBMHX(760:770,:));

corrInd=medianHXend>-100;

data.LDR(data.FLAG>1)=nan;

dbm=data.DBMHX(:,corrInd==1);
ldr=data.LDR(:,corrInd==1);

%% Find filter indices
ldrMask=~isnan(ldr);

ldrMaskOut=zeros(size(ldrMask));

for ii=1:size(ldrMask,2)
    rayIn=ldrMask(:,ii);    
    ldrMaskOut(:,ii)=bwareaopen(rayIn,10);
end

ldrMaskOut=bwareaopen(ldrMaskOut,20);

filterLDR=(ldrMask==1 & ldrMaskOut==0);
%% Plot
ylimits=[-0.5 15];

close all

figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200]);

ax1=subplot(3,1,1);
fig1=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
%title('Reflectivity (all data)')

ax1=subplot(3,1,2);
fig1=surf(data.time,data.asl./1000,double(filterLDR),'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
%title('Reflectivity (all data)')
% 
% linkaxes([ax1,ax2,ax3],'xy');

% formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
% print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_echoMask'],...
%    '-dpng','-r0');
% 
