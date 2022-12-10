% SST on flight path

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

cfDir=HCRdir(project,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

savedir=[cfDir(1:end-5),'cloudProps/'];

%% Initiate output

sstAll=[];
lonAll=[];
latAll=[];

%% Loop through flights

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    fileList=makeFileList(cfDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];
    data.SST=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    sstAll=cat(2,sstAll,data.SST);
    lonAll=cat(2,lonAll,data.longitude);
    latAll=cat(2,latAll,data.latitude);
end

%% Save properties

disp('Saving output ...');

save([savedir,project,'_ssts.mat'],'sstAll','lonAll','latAll');
