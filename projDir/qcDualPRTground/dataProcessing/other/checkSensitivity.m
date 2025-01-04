% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'sensitivity/'];

%% Run processing
infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];
caseList = table2array(readtable(infile));

%% Run processing

longSNR=[];
shortSNR=[];
longDBZ=[];
shortDBZ=[];

% Go through flights
for ii=1:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];
    data.SNRVC_short=[];
    data.DBZ_short=[];
    data.FLAG_short=[];
    data.SNRVC_long=[];
    data.DBZ_long=[];
    data.FLAG_long=[];

    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.SNRVC_short(data.FLAG_short~=1)=nan;
    data.SNRVC_long(data.FLAG_long~=1)=nan;
    data.DBZ_short(data.FLAG_short~=1)=nan;
    data.DBZ_long(data.FLAG_long~=1)=nan;

    longNan=all(isnan(data.SNRVC_long),1);

    longSNR=cat(2,longSNR,data.SNRVC_long(:,longNan==0));
    shortSNR=cat(2,shortSNR,data.SNRVC_short(:,longNan==0));
    longDBZ=cat(2,longDBZ,data.DBZ_long(:,longNan==0));
    shortDBZ=cat(2,shortDBZ,data.DBZ_short(:,longNan==0));
end

range1=data.range(:,1);
save([figdir,'SNRs.mat'],'longSNR','shortSNR','longDBZ','shortDBZ','range1','-v7.3');