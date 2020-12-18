% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

startTime=datetime(2018,1,19,5,30,0);
endTime=datetime(2018,1,19,6,45,0);

%% Load data

disp('Loading data ...');

data=[];

data.MELTING_LAYER=[];
data.ICING_LEVEL=[];

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

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);
%% Plot

timeMat=repmat(data.time,size(data.MELTING_LAYER,1),1);


close all

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end

% Resample for plotting
newFindMelt=data.MELTING_LAYER(:,newInds);
newASL=data.asl(:,newInds);
newTime=data.time(newInds);

%% Plot
close all

%fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);
wi=5;
hi=3;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');
colmap=jet;
colormap(flipud(colmap));

ax2=subplot(1,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
ax2.Colormap=([1 0 1;1 1 0]);
view(2);
scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');

scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');
scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');

plot(data.time,data.ICING_LEVEL./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
ax = gca;
ax.SortMethod = 'childorder';
ylim([0 4]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Melting layer data product')
grid on
ax2.Position=[0.085 0.157 0.87 0.76];

print([figdir,'icinglevel'],'-dpng','-r0');

