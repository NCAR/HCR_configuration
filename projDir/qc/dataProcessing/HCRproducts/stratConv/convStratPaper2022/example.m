% Plot HCR convStrat from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';
whichModel='era5';

% startTime=datetime(2019,9,25,15,12,0);
% endTime=datetime(2019,9,25,15,34,0);

startTime=datetime(2019,9,27,12,32,0);
endTime=datetime(2019,9,27,12,50,50);

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    0.7,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    1,0,0];

colmapBasic=[0,0,1;
    0.32,0.78,0.59;
    0.7,0,0];


%% Get data

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ_MASKED=[];
data.VEL_MASKED=[];
data.CONVECTIVITY=[];
data.PARTITION_2D=[];
data.PARTITION_1D=[];

dataVars=fieldnames(data);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Plot in hourly increments

disp('Plotting ...');

startPlot=startTime;
endPlot=endTime;

timePlot=data.time;
dbzPlot=data.DBZ_MASKED;

aslPlot=data.asl;
velPlot=data.VEL_MASKED;

sc1D=data.PARTITION_1D;
sc1D(sc1D==14)=1;
sc1D(sc1D==16)=2;
sc1D(sc1D==18)=3;
sc1D(sc1D==25)=4;
sc1D(sc1D==30)=5;
sc1D(sc1D==32)=6;
sc1D(sc1D==34)=7;
sc1D(sc1D==36)=8;
sc1D(sc1D==38)=9;

convPlot=data.CONVECTIVITY;

time1D=timePlot(~isnan(sc1D));
sc1D=sc1D(~isnan(sc1D));
col1D=colmapSC(sc1D,:);

stratConvPlot=data.PARTITION_2D;
stratConvPlot(stratConvPlot==14)=1;
stratConvPlot(stratConvPlot==16)=2;
stratConvPlot(stratConvPlot==18)=3;
stratConvPlot(stratConvPlot==25)=4;
stratConvPlot(stratConvPlot==30)=5;
stratConvPlot(stratConvPlot==32)=6;
stratConvPlot(stratConvPlot==34)=7;
stratConvPlot(stratConvPlot==36)=8;
stratConvPlot(stratConvPlot==38)=9;

f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible','on');

s1=subplot(4,1,1);

colormap jet

hold on
surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
colorbar
grid on
title('Reflectivity (dBZ)')
s1pos=s1.Position;
s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];

s2=subplot(4,1,2);

colormap jet

hold on
surf(timePlot,aslPlot./1000,velPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-5 5]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
colorbar
grid on
title('Velocity (m s^{-1})')
s2pos=s2.Position;
s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

s3=subplot(4,1,3);

colormap jet

hold on
surf(timePlot,aslPlot./1000,convPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
colorbar
grid on
title('Convectivity');
s3pos=s3.Position;
s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

s5=subplot(30,1,30);

hold on
scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
%set(gca,'clim',[0,1]);
set(gca,'YTickLabel',[]);
s5.Colormap=colmapSC;
xlim([timePlot(1),timePlot(end)]);
s5pos=s5.Position;
s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];

s4=subplot(4,1,4);

hold on
surf(timePlot,aslPlot./1000,stratConvPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 10]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
s4.Colormap=colmapSC;
caxis([0.5 9.5]);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);
grid on
title('Stratiform/convective partitioning')
s4pos=s4.Position;
s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'algorithm.png'],'-dpng','-r0')
