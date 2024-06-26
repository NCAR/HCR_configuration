% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='';
qcVersion='';

startTime=datetime(2024,6,24,22,31,0);
endTime=datetime(2024,6,24,22,32,0);

dataDirTS=HCRdir(project,quality,qcVersion,freqData);

%% Get data

disp('Reading data ...');

fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

        data=[];
            data.IVc=[];
            data.QVc=[];
            data=read_TsArchive_iwrf_bulk(fileListTS{1},data);




