% Plot HCR variables

clear all
close all

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.0';

startTime=datetime(2021,6,1,20,46,0);
endTime=datetime(2021,6,1,20,56,0);

ylimUpper=11;

saveFig=0;
figdir=['.'];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

%% Get data

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBMVC=[];
data.DBMHX=[];
data.SNRHX=[];
data.LDR=[];

data=read_HCR(fileList,data,startTime,endTime);

[fig,s]=do_plotHCR(data,ylimUpper);

s.DBMVC.CLim=[-110,-90];
s.DBMHX.CLim=[-102.5,-101.5];
s.SNRHX.CLim=[-15,-10];
s.LDR.CLim=[-25,0];