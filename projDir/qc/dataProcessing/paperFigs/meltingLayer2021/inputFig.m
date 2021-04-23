% Analyze HCR clouds
clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.1';
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));
%figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];
figdir='/scr/sci/romatsch/paperFigs/meltLayer/';

%indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];
indir=HCRdir(project,quality,qcVersion,freqData);

startTime=datetime(2018,2,5,0,5,0);
endTime=datetime(2018,2,5,0,37,0);

%% Load data

disp('Loading data ...');

data=[];

data.TEMP=[];
data.VEL_CORR=[];
data.LDR=[];
data.FLAG=[];
data.MELTING_LAYER=[];

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
twentyoneInds=find(data.MELTING_LAYER==21);

%% Plot

timeMat=repmat(data.time,size(data.LDR,1),1);
ldrMasked=data.LDR;
ldrMasked(data.FLAG>1)=nan;
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
newTEMP=data.TEMP(:,newInds);
newLDR=ldrMasked(:,newInds);
newVEL=velMasked(:,newInds);
newASL=data.asl(:,newInds);
newTime=data.time(newInds);

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

ylimits=[0 4];

ax1=subplot(3,1,1);

hold on;
sub1=surf(newTime,newASL./1000,newTEMP,'edgecolor','none');
view(2);
colormap(jet)
caxis([-10 10]);
hcb1=colorbar;
l1=scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');

scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');

ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title(['(a) TEMP (',char(176),'C)'])
grid on
ax1.Position=[0.08 0.71 0.8 0.26];
hcb1.Position=[0.895 0.71 0.04 0.26];
legend(l1,{'Zero deg'},'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2=subplot(3,1,2);
hold on;
sub3=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
view(2);
caxis([-25 -5]);
hcb2=colorbar;
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(b) LDR (dB)')
grid on
ax2.Position=[0.08 0.385 0.8 0.26];
hcb2.Position=[0.895 0.385 0.04 0.26];

ax3=subplot(3,1,3);
ax3.Colormap=jet;
hold on;
sub3=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
caxis([-4 4]);
hcb3=colorbar('XTick',-4:2:4);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('(c) VEL (m s^{-1})')
grid on
ax3.Position=[0.08 0.06 0.8 0.26];
hcb3.Position=[0.895 0.06 0.04 0.26];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'inputFig'],'-dpng','-r0');
