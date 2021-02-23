% Analyze HCR clouds
clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

startTime=datetime(2019,8,7,17,26,0);
endTime=datetime(2019,8,7,17,32,0);

%% Load data

disp('Loading data ...');

data=[];

data.MELTING_LAYER=[];
data.ICING_LEVEL=[];
data.VEL_CORR=[];
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

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);

%% Prepare

timeMat=repmat(data.time,size(data.VEL_CORR,1),1);
velMasked=data.VEL_CORR;
velMasked(data.FLAG>1)=nan;

close all

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end

% Resample for plotting
newVEL=velMasked(:,newInds);
newASL=data.asl(:,newInds);
newFindMelt=data.MELTING_LAYER(:,newInds);
newTime=data.time(newInds);

%% Plot

close all

wi=5;
hi=6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');
colormap(jet);

ylimits=[2 8];

ax1=subplot(2,1,1);

hold on;
sub1=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(a) Radial velocity (m s^{-1})')
grid on
caxis([-5 5]);
hcb=colorbar('XTick',-4:2:4);
ax1.Position=[0.11 0.577 0.77 0.38];
hcb.Position=[0.895 0.58 0.04 0.34];

ax2=subplot(2,1,2);
hold on;
sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
ax2.Colormap=([1 0 1;1 1 0]);
view(2);
scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,7,'k','filled','MarkerEdgeColor','k');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l1=scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,7,'k','filled','MarkerEdgeColor','k');
l2=scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
l3=scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
l4=scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l5=plot(data.time,data.ICING_LEVEL./1000,'linewidth',1.2,'color',[0.8 0.8 0.8]);
ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(b) Melting layer data product')
grid on
legend([l1 l2 l3 l4 l5],{'Zero deg','Estimates','Interpolations','Detections','Icing level'},...
    'Location','northeast');
ax2.Position=[0.11 0.08 0.77 0.38];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer'],'-dpng','-r0');
