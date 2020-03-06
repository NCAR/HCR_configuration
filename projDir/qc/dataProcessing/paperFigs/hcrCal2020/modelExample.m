% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2018,1,24,0,45,0);
endTime=datetime(2018,1,24,1,30,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

directories.dataDir=HCRdir(project,quality,freqData);

if strcmp(whichModel,'era5')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/era5interp/',project,'/',freqData,'/'];
elseif strcmp(whichModel,'ecmwf')
    directories.modeldir=['/scr/sci/romatsch/data/reanalysis/ecmwf/forecastInterp/',project,'/',freqData,'/'];
end

%% Load data

data.DBZ=[];
data.RH=[];
data.TEMP=[];

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

upInd=find(data.elevation>0);

%% Find zero degree altitude
tempVec=1:1:size(data.TEMP,1);
tempMat=repmat(tempVec,size(data.TEMP,2),1)';

tempMat(data.TEMP>0)=nan;
tempMat(isnan(data.TEMP))=nan;

[maxs,rowInds] = nanmax(tempMat);
colInds=1:1:size(data.TEMP,2);

linearInd = sub2ind(size(data.TEMP), rowInds, colInds);

[mins,rowIndsUp] = nanmin(tempMat);
colIndsUp=1:1:size(data.TEMP,2);

linearIndUp = sub2ind(size(data.TEMP), rowIndsUp, colIndsUp);

linearInd(upInd)=linearIndUp(upInd);
rowInds(upInd)=rowIndsUp(upInd);

checkTemps=data.TEMP(linearInd);
outInds=zeros(size(checkTemps));
outInds(checkTemps>0 | checkTemps<-2)=1;
outInds(isnan(checkTemps))=1;
linearInd(outInds==1)=[];
zTempsNeg=data.asl(linearInd);

timeMat=repmat(data.time,size(data.TEMP,1),1);

%% Plot
ylimits=[-0.1 6];

close all

colHum=jet;
colHum=flipud(colHum(30:end,:));

wi=10;
hi=6.5;

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
ax1=subplot(2,1,1);
hold on;
outerpos1 = ax1.Position;
ax1.Position = [outerpos1(1)-0.07 outerpos1(2)+0.0 outerpos1(3)+0.14 outerpos1(4)+0.02];

fig1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
fig1=colMapDBZ(fig1);
plot(timeMat(linearInd),zTempsNeg./1000,'c','linewidth',2);
ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(a)  Reflectivity and ERA5 freezing level')

ax2=subplot(2,1,2);
hold on;
outerpos1 = ax2.Position;
ax2.Position = [outerpos1(1)-0.07 outerpos1(2)-0.02 outerpos1(3)+0.14 outerpos1(4)+0.02];
fig2=surf(data.time,data.asl./1000,data.RH,'edgecolor','none');
view(2);
caxis([0 100]);
colormap(ax2,colHum);
hcb=colorbar;
set(get(hcb,'Title'),'String','%');
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(b)                      ERA5 relative humidity')

print([figdir,'modelExample.png'],'-dpng','-r0');