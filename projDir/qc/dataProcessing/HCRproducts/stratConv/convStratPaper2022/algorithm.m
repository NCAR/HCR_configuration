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
data.ECHO_TYPE_2D=[];
data.ECHO_TYPE_1D=[];

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

sc1D=data.ECHO_TYPE_1D;
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

stratConvPlot=data.ECHO_TYPE_2D;
stratConvPlot(stratConvPlot==14)=1;
stratConvPlot(stratConvPlot==16)=2;
stratConvPlot(stratConvPlot==18)=3;
stratConvPlot(stratConvPlot==25)=4;
stratConvPlot(stratConvPlot==30)=5;
stratConvPlot(stratConvPlot==32)=6;
stratConvPlot(stratConvPlot==34)=7;
stratConvPlot(stratConvPlot==36)=8;
stratConvPlot(stratConvPlot==38)=9;

stratConvPlotBasic=nan(size(stratConvPlot));
stratConvPlotBasic(stratConvPlot<4)=1;
stratConvPlotBasic(stratConvPlot==4)=2;
stratConvPlotBasic(stratConvPlot>4)=3;


%% Plot

close all

wi=9;
hi=10;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(6,1,1);

colormap jet

hold on
surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
set(gca,'XTickLabel',[]);
cb1=colorbar;
plot(data.time,data.altitude./1000,'-k','LineWidth',2);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

s2=subplot(6,1,2);

colormap jet

hold on
surf(timePlot,aslPlot./1000,velPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-5 5]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
set(gca,'XTickLabel',[]);
cb2=colorbar;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(b) Velocity (m s^{-1})','FontSize',11,'FontWeight','bold');

s3=subplot(6,1,3);

colormap jet

hold on
surf(timePlot,aslPlot./1000,convPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
set(gca,'XTickLabel',[]);
cb3=colorbar;
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(c) Convectivity','FontSize',11,'FontWeight','bold');

s4=subplot(6,1,4);

hold on
surf(timePlot,aslPlot./1000,stratConvPlotBasic,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
s4.Colormap=colmapBasic;
caxis([0.5 3.5]);
cb4=colorbar;
cb4.Ticks=1:3;
cb4.TickLabels={'Stratiform','Mixed','Convective'};
set(gca,'XTickLabel',[]);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(d) Basic echo type','FontSize',11,'FontWeight','bold');

s6=subplot(35,1,35);

hold on
scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
%set(gca,'clim',[0,1]);
set(gca,'YTickLabel',[]);
s6.Colormap=colmapSC;
xlim([timePlot(1),timePlot(end)]);
grid on
box on

s5=subplot(6,1,5);

hold on
surf(timePlot,aslPlot./1000,stratConvPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
ylim([0 ylimUpper]);
xlim([timePlot(1),timePlot(end)]);
s5.Colormap=colmapSC;
caxis([0.5 9.5]);
cb5=colorbar;
cb5.Ticks=1:9;
cb5.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);
grid on
box on
text(datetime(2019,9,27,12,32,10),14,'(e) Advanced echo type','FontSize',11,'FontWeight','bold');

s1.Position=[0.049 0.805 0.82 0.17];
s2.Position=[0.049 0.62 0.82 0.17];
s3.Position=[0.049 0.435 0.82 0.17];
s4.Position=[0.049 0.25 0.82 0.17];
s5.Position=[0.049 0.065 0.82 0.17];
s6.Position=[0.049 0.04 0.82 0.016];

cb1.Position=[0.875,0.805,0.02,0.17];
cb2.Position=[0.875,0.62,0.02,0.17];
cb3.Position=[0.875,0.435,0.02,0.17];
cb4.Position=[0.875,0.25,0.02,0.17];
cb5.Position=[0.875,0.065,0.02,0.17];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'algorithm.png'],'-dpng','-r0')
