% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.1';

startTime=datetime(2021,5,29,18,50,0);
endTime=datetime(2021,5,29,19,32,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=12;

saveFig=1;
if saveFig
    outname='spiculeDBZ'
    %figdir=['.'];
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

[fig,s]=do_plotHCR(data,ylimUpper);

%s.DBMHX.CLim=[-103,-90];

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end