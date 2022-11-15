% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

startTime=datetime(2019,8,16,14,52,0);
endTime=datetime(2019,8,16,15,8,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=14;

saveFig=1;
if saveFig
    outname=['echoType_',datestr(startTime,'yyyymmdd_HHMMSS')];
    %figdir=['/scr/sci/romatsch/HCR/examplePlots/'];
    figdir=[indir(1:end-5),'examplePlots/'];
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
end

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.ECHO_TYPE_2D=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

fig=figure('Position',[200 500 1200 600],'DefaultAxesFontSize',14);

classSubPlot=data.ECHO_TYPE_2D;

classSubPlot(classSubPlot==14)=1;
classSubPlot(classSubPlot==16)=2;
classSubPlot(classSubPlot==18)=3;
classSubPlot(classSubPlot==25)=4;
classSubPlot(classSubPlot==30)=5;
classSubPlot(classSubPlot==32)=6;
classSubPlot(classSubPlot==34)=7;
classSubPlot(classSubPlot==36)=8;
classSubPlot(classSubPlot==38)=9;


surf(data.time,data.asl./1000,classSubPlot,'EdgeColor','none');
view(2)
caxis([1,12]);

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    0.7,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    1,0,0];

fig.Colormap=colmapSC;
caxis([0.5 9.5]);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};

title('Echo type')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end