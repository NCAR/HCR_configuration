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
    outname=['flag_',datestr(startTime,'yyyymmdd_HHMMSS')];
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

data.FLAG=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

fig=figure('Position',[200 500 1200 600],'DefaultAxesFontSize',14);

surf(data.time,data.asl./1000,data.FLAG,'EdgeColor','none');
view(2)
caxis([1,12]);

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

colMask=[0.4,0.8,1;
        0,0,0;
        0.5,0,0.5;
        0,1,0;
        0.2,0.6,0.2;
        1,0,0;
        0,0,0.6;
        0.7065,0.4415,0.2812;
        0.5,0.5,0.5;
        0.9290,0.8940,0.1250;
        1,0,1;
        1,0.6,0];

ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
        'Bang (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
        'NS cal (10)';'Ant. trans. (11)';'Missing (12)'};

fig.Colormap=colMask;
        hcb=colorbar;
        hcb.Ticks=[1.5:0.91:12.5];
        hcb.TickLabels=ytickLabels;

title('Flag')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end