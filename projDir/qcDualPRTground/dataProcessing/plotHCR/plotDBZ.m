% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset, otrec
quality='qc0'; %field, qc1, or qc2
freqData='100hz_long';
qcVersion='';

startTime=datetime(2024,6,24,22,31,0);
endTime=datetime(2024,6,24,22,39,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=12;

saveFig=1;
if saveFig
    outname=['dbz_',datestr(startTime,'yyyymmdd_HHMMSS')];
    %figdir=['/scr/sci/romatsch/HCR/examplePlots/'];
    figdir=[indir(1:end-14),'examplePlots/'];
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
end

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ_long=[];

data=read_HCR(fileList,data,startTime,endTime);


%% Plot
pix=1;

close all

disp('Plotting ...');

fig=figure('Position',[200 500 1200 600],'DefaultAxesFontSize',14);

s1=surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBZ_long(:,1:pix:length(data.time)),'EdgeColor','none');
view(2)
clim([-40 25]);
colM=jet;
colormap(fig,colM);
colorbar

xlim([data.time(1),data.time(end)]);
ylim([0 ylimUpper]);

ylabel('Range (km)')

grid on
box on

title('Reflectivity (dBZ)')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end