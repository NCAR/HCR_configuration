% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

startTime=datetime(2018,1,23,0,28,0);
endTime=datetime(2018,1,23,0,38,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=4;

saveFig=1;
if saveFig
    outname='socratesCONVECTIVITY';
    figdir=['/scr/sci/romatsch/HCR/examplePlots/'];
    %figdir=[indir(1:end-5),'examplePlots/'];
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
end

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.CONVECTIVITY=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

[fig,s]=do_plotHCR(data,ylimUpper);

s.VEL_MASKED.CLim=[-7,7];

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end