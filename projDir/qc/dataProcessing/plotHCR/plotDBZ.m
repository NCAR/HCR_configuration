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
    outname=['dbz_',datestr(startTime,'yyyymmdd_HHMMSS')];
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

data.DBZ_MASKED=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

fig=figure('Position',[200 500 1200 600],'DefaultAxesFontSize',14);

s1=surf(data.time,data.asl./1000,data.DBZ_MASKED,'EdgeColor','none');
view(2)
colorbar
s1=colMapDBZ(s1);

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

title('Reflectivity (dBZ)')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end