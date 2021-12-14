% Plot HCR variables

clear all
close all

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.0';

startTime=datetime(2021,6,17,19,00,0);
endTime=datetime(2021,6,17,19,40,0);

ylimUpper=14;

saveFig=0;
figdir=['.'];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBMHX=[];

data=read_HCR(fileList,data,startTime,endTime);

disp('Plotting ...');

[fig,s]=do_plotHCR(data,ylimUpper);

s.DBMHX.CLim=[-103,-90];
