% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

%% Load data1

disp('Loading data 1 ...');

startTime=datetime(2018,1,19,5,30,0);
endTime=datetime(2018,1,19,6,45,0);

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

timeMat=repmat(data.time,size(data.MELTING_LAYER,1),1);

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

%% Load data2

disp('Loading data 2 ...');

startTime=datetime(2018,2,4,3,30,0);
endTime=datetime(2018,2,4,4,50,0);

data2=[];

data2.MELTING_LAYER=[];
data2.ICING_LEVEL=[];

dataVars=fieldnames(data2);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data2=read_HCR(fileList,data2,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data2,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

elevenInds2=find(data2.MELTING_LAYER==11);
twelveInds2=find(data2.MELTING_LAYER==12);
thirteenInds2=find(data2.MELTING_LAYER==13);
fourteenInds2=find(data2.MELTING_LAYER==14);

twentyoneInds2=find(data2.MELTING_LAYER==21);
twentytwoInds2=find(data2.MELTING_LAYER==22);
twentythreeInds2=find(data2.MELTING_LAYER==23);
twentyfourInds2=find(data2.MELTING_LAYER==24);

timeMat2=repmat(data2.time,size(data2.MELTING_LAYER,1),1);

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data2.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data2.time);
else
    newInds=1:100:length(data2.time);
end

% Resample for plotting
newFindMelt2=data2.MELTING_LAYER(:,newInds);
newASL2=data2.asl(:,newInds);
newTime2=data2.time(newInds);

%% Load data3

disp('Loading data 3 ...');

startTime=datetime(2018,2,20,4,10,0);
endTime=datetime(2018,2,20,4,12,0);

data3=[];

data3.MELTING_LAYER=[];
data3.ICING_LEVEL=[];
data3.LDR=[];
data3.FLAG=[];

dataVars=fieldnames(data3);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data3=read_HCR(fileList,data3,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data3,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

elevenInds3=find(data3.MELTING_LAYER==11);
twelveInds3=find(data3.MELTING_LAYER==12);
thirteenInds3=find(data3.MELTING_LAYER==13);
fourteenInds3=find(data3.MELTING_LAYER==14);

twentyoneInds3=find(data3.MELTING_LAYER==21);
twentytwoInds3=find(data3.MELTING_LAYER==22);
twentythreeInds3=find(data3.MELTING_LAYER==23);
twentyfourInds3=find(data3.MELTING_LAYER==24);

timeMat3=repmat(data3.time,size(data3.MELTING_LAYER,1),1);

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data3.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data3.time);
else
    newInds=1:100:length(data3.time);
end

data3.LDR(data3.FLAG>1)=nan;

% Resample for plotting
newFindMelt3=data3.MELTING_LAYER(:,newInds);
newASL3=data3.asl(:,newInds);
newLDR3=data3.LDR(:,newInds);
newTime3=data3.time(newInds);

%% Plot
close all

wi=10;
hi=6.5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

ax1=subplot(2,2,1);
hold on;
surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
ax1.Colormap=([1 0 1;1 1 0]);
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
ylim([0 4]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(a) Melting layer data product')
grid on
ax1.Position=[0.048 0.57 0.43 0.38];

ax3=subplot(2,2,3);
hold on;
surf(newTime2,newASL2./1000,newFindMelt2,'edgecolor','none');
ax3.Colormap=([1 0 1;1 1 0]);
view(2);
scatter(timeMat2(elevenInds2),data2.asl(elevenInds2)./1000,7,'k','filled','MarkerEdgeColor','k');
scatter(timeMat2(fourteenInds2),data2.asl(fourteenInds2)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(timeMat2(thirteenInds2),data2.asl(thirteenInds2)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
scatter(timeMat2(twelveInds2),data2.asl(twelveInds2)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l1=scatter(timeMat2(twentyoneInds2),data2.asl(twentyoneInds2)./1000,7,'k','filled','MarkerEdgeColor','k');
l2=scatter(timeMat2(twentyfourInds2),data2.asl(twentyfourInds2)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
l3=scatter(timeMat2(twentythreeInds2),data2.asl(twentythreeInds2)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
l4=scatter(timeMat2(twentytwoInds2),data2.asl(twentytwoInds2)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l5=plot(data2.time,data2.ICING_LEVEL./1000,'linewidth',1.2,'color',[0.8 0.8 0.8]);

ax = gca;
ax.SortMethod = 'childorder';
ylim([0 2]);
ylabel('Altitude (km)');
xlim([data2.time(1),data2.time(end)]);
title('(b) Melting layer data product')
grid on
yticks(0:1:2)
legend([l1 l2 l3 l4 l5],{'Zero deg','Estimates','Interpolations','Detections','Icing level'},...
    'Location','northwest');
ax3.Position=[0.048 0.07 0.43 0.38];

ax2=subplot(2,2,2);
hold on;
sub1=surf(newTime3,newASL3./1000,newLDR3,'edgecolor','none');
ax2.Colormap=(jet);
hcb2=colorbar;
view(2);
ylim([0 2]);
yticks(0:1:2)
caxis([-35 -5]);
ylabel('Altitude (km)');
xlim([data3.time(1),data3.time(end)]);
title('(c) LDR (dB)')
grid on
ax2.Position=[0.535 0.57 0.41 0.38];
hcb2.Position=[0.95 0.57 0.02 0.38];

ax4=subplot(2,2,4);
hold on;
sub1=surf(newTime3,newASL3./1000,newFindMelt3,'edgecolor','none');
ax4.Colormap=([1 0 1;1 1 0]);
view(2);
scatter(timeMat3(elevenInds3),data3.asl(elevenInds3)./1000,7,'k','filled','MarkerEdgeColor','k');
scatter(timeMat3(fourteenInds3),data3.asl(fourteenInds3)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(timeMat3(thirteenInds3),data3.asl(thirteenInds3)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
scatter(timeMat3(twelveInds3),data3.asl(twelveInds3)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l1=scatter(timeMat3(twentyoneInds3),data3.asl(twentyoneInds3)./1000,7,'k','filled','MarkerEdgeColor','k');
l2=scatter(timeMat3(twentyfourInds3),data3.asl(twentyfourInds3)./1000,7,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
l3=scatter(timeMat3(twentythreeInds3),data3.asl(twentythreeInds3)./1000,7,'c','filled','filled','MarkerEdgeColor','c');
l4=scatter(timeMat3(twentytwoInds3),data3.asl(twentytwoInds3)./1000,7,'b','filled','filled','MarkerEdgeColor','b');

l5=plot(data3.time,data3.ICING_LEVEL./1000,'linewidth',1.2,'color',[0.8 0.8 0.8]);

ax = gca;
ax.SortMethod = 'childorder';
ylim([0 2]);
yticks(0:1:2)
ylabel('Altitude (km)');
xlim([data3.time(1),data3.time(end)]);
title('(d) Melting layer data product')
grid on
ax4.Position=[0.535 0.07 0.41 0.38];

print([figdir,'challenges'],'-dpng','-r0');

