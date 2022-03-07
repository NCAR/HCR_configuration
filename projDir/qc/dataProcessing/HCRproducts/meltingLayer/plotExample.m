% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='combined'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
if strcmp(project,'otrec')
    ylimits=[-0.3 14];
elseif strcmp(project,'socrates')
    ylimits=[-0.2 6];
elseif strcmp(project,'cset')
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/meltLayerHCR/',project,'/cases/'];
%figdir=['/home/romatsch/plots/HCR/meltingLayer/selected/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/meltLayer_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);
%indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];


startTime=datetime(2018,1,22,22,47,0);
endTime=datetime(2018,1,22,23,1,0);

%% Load data

disp('Loading data ...');

data=[];

data.HCR_DBZ=[];
data.HCR_MELTING_LAYER=[];
data.HCR_ECHO_TYPE_2D=[];
data.PID=[];

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

data.DBZ=data.HCR_DBZ;
data.MELTING_LAYER=data.HCR_MELTING_LAYER;

%% Find melting layer

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);

%% Plot

timeMat=repmat(data.time,size(data.DBZ,1),1);

close all

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    0.7,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    1,0,0];

stratConvPlot=data.HCR_ECHO_TYPE_2D;
stratConvPlot(stratConvPlot==14)=1;
stratConvPlot(stratConvPlot==16)=2;
stratConvPlot(stratConvPlot==18)=3;
stratConvPlot(stratConvPlot==25)=4;
stratConvPlot(stratConvPlot==30)=5;
stratConvPlot(stratConvPlot==32)=6;
stratConvPlot(stratConvPlot==34)=7;
stratConvPlot(stratConvPlot==36)=8;
stratConvPlot(stratConvPlot==38)=9;

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];

units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
    'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,800]);

ax1=subplot(3,1,1);
hold on;
sub1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
sub1=colMapDBZ(sub1);
scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');

scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');
scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');
ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity and melting layer')
grid on
set(gca,'xticklabel',[])

ax2=subplot(3,1,2);
hold on;
sub1=surf(data.time,data.asl./1000,stratConvPlot,'edgecolor','none');
caxis([0.5 9.5]);
cb2=colorbar;
cb2.Ticks=1:9;
cb2.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
ax2.Colormap=colmapSC;
view(2);

ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Convective/stratiform echo type')
grid on
set(gca,'xticklabel',[])

%%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3=subplot(3,1,3);
hold on;
sub3=surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
view(2);
colormap(ax3,cscale_hcr);
cb=colorbar;
cb.Ticks=1:11;
cb.TickLabels=units_str_hcr;
caxis([.5 11.5]);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Hydrometeor type')
grid on

ax1.Position=[0.035 0.7 0.79 0.27];
ax2.Position=[0.035 0.38 0.79 0.27];
ax3.Position=[0.035 0.06 0.79 0.27];

linkaxes([ax1 ax2 ax3],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
