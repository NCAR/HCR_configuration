% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

%% Load data1

disp('Loading data 1 ...');

startTime=datetime(2019,9,21,12,25,0);
endTime=datetime(2019,9,21,12,29,0);

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

timeMat=repmat(data.time,size(data.MELTING_LAYER,1),1);

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end
data.VEL_CORR(data.FLAG>1)=nan;

% Resample for plotting
newFindMelt=data.MELTING_LAYER(:,newInds);
newASL=data.asl(:,newInds);
newVEL=data.VEL_CORR(:,newInds);
newTime=data.time(newInds);

%% Load data2

disp('Loading data 2 ...');

startTime=datetime(2019,8,16,15,34,0);
endTime=datetime(2019,8,16,15,42,0);

data2=[];

data2.MELTING_LAYER=[];
data2.ICING_LEVEL=[];
data2.VEL_CORR=[];
data2.FLAG=[];

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

data2.VEL_CORR(data2.FLAG>1)=nan;

% Resample for plotting
newFindMelt2=data2.MELTING_LAYER(:,newInds);
newASL2=data2.asl(:,newInds);
newVEL2=data2.VEL_CORR(:,newInds);
newTime2=data2.time(newInds);

%% Load data3

disp('Loading data 3 ...');

startTime=datetime(2019,8,25,16,17,0);
endTime=datetime(2019,8,25,16,19,0);

data3=[];

data3.MELTING_LAYER=[];
data3.ICING_LEVEL=[];
data3.VEL_CORR=[];
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

data3.VEL_CORR(data3.FLAG>1)=nan;

% Resample for plotting
newFindMelt3=data3.MELTING_LAYER(:,newInds);
newASL3=data3.asl(:,newInds);
newVEL3=data3.VEL_CORR(:,newInds);
newTime3=data3.time(newInds);

%% Plot
close all

wi=5;
hi=8;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

ax1=subplot(3,1,1);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
ax1.Colormap=(jet);
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
ylim([4 5]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(a) VEL (m s^{-1})')
grid on
caxis([-6 6]);
hcb1=colorbar('XTick',-6:2:6);
ax1.Position=[0.11 0.71 0.78 0.26];
hcb1.Position=[0.9 0.71 0.04 0.26];

ax2=subplot(3,1,2);
hold on;
surf(newTime2,newASL2./1000,newVEL2,'edgecolor','none');
ax2.Colormap=(jet);
view(2);
l1=scatter(timeMat2(elevenInds2),data2.asl(elevenInds2)./1000,10,'k','filled');
l2=scatter(timeMat2(fourteenInds2),data2.asl(fourteenInds2)./1000,10,'g','filled');
l3=scatter(timeMat2(thirteenInds2),data2.asl(thirteenInds2)./1000,10,'c','filled');
l4=scatter(timeMat2(twelveInds2),data2.asl(twelveInds2)./1000,10,'b','filled');

scatter(timeMat2(twentyoneInds2),data2.asl(twentyoneInds2)./1000,10,'k','filled');
scatter(timeMat2(twentyfourInds2),data2.asl(twentyfourInds2)./1000,10,'g','filled');
scatter(timeMat2(twentythreeInds2),data2.asl(twentythreeInds2)./1000,10,'c','filled');
scatter(timeMat2(twentytwoInds2),data2.asl(twentytwoInds2)./1000,10,'b','filled');

l5=plot(data2.time,data2.ICING_LEVEL./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
ax = gca;
ax.SortMethod = 'childorder';
ylim([4 5]);
ylabel('Altitude (km)');
xlim([data2.time(1),data2.time(end)]);
title('(b) VEL (m s^{-1})')
grid on
caxis([-6 6]);
hcb2=colorbar('XTick',-6:2:6);
ax2.Position=[0.11 0.385 0.78 0.26];
hcb2.Position=[0.9 0.385 0.04 0.26];

ax3=subplot(3,1,3);
hold on;
sub1=surf(newTime3,newASL3./1000,newVEL3,'edgecolor','none');
ax3.Colormap=(jet);
view(2);
scatter(timeMat3(elevenInds3),data3.asl(elevenInds3)./1000,10,'k','filled');
scatter(timeMat3(fourteenInds3),data3.asl(fourteenInds3)./1000,10,'g','filled');
scatter(timeMat3(thirteenInds3),data3.asl(thirteenInds3)./1000,10,'c','filled');
scatter(timeMat3(twelveInds3),data3.asl(twelveInds3)./1000,10,'b','filled');

scatter(timeMat3(twentyoneInds3),data3.asl(twentyoneInds3)./1000,10,'k','filled');
scatter(timeMat3(twentyfourInds3),data3.asl(twentyfourInds3)./1000,10,'g','filled');
scatter(timeMat3(twentythreeInds3),data3.asl(twentythreeInds3)./1000,10,'c','filled');
scatter(timeMat3(twentytwoInds3),data3.asl(twentytwoInds3)./1000,10,'b','filled');

plot(data3.time,data3.ICING_LEVEL./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
ax = gca;
ax.SortMethod = 'childorder';
ylim([4 5]);
%yticks(0:1:2)
hcb1=colorbar('XTick',-4:2:4);
ylabel('Altitude (km)');
xlim([data3.time(1),data3.time(end)]);
title('(c) VEL (m s^{-1})')
grid on
caxis([-6 6]);
hcb3=colorbar('XTick',-6:2:6);
ax3.Position=[0.11 0.06 0.78 0.26];
hcb3.Position=[0.9 0.06 0.04 0.26];

print([figdir,'variability'],'-dpng','-r0');

