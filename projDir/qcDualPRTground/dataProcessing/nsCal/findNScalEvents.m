% Find nsCal events
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; % socrates, cset, aristo, otrec
quality='qc0'; % field, qc0, qc1, qc2
qcVersion='v1.0';
freq='10hz_combined'; %10hz_combined, 50hz_longPulse, or 50hz_shortPulse

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freq);

%% Run processing

startEndAll=[];

% Go through flights
for ii=1:size(caseList,1)
    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    data=[];
    data.DBMVC_short=[];
        
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
       
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
        
    %% Find nsCal
    firstGate=data.DBMVC_short(1,:);
    nscalInds=find(firstGate>-96 & firstGate<-87);
    nscalMask=zeros(size(firstGate));
    nscalMask(nscalInds)=1;
    
    diffMask=diff(nscalMask);
    
    startNSinds=find(diffMask==1);
    endNSinds=find(diffMask==-1);
    
    startEndDiff=endNSinds-startNSinds;
    
    startNSinds(startEndDiff<60)=[];
    endNSinds(startEndDiff<60)=[];
    
    startNStime=data.time(startNSinds);
    endNStime=data.time(endNSinds);
    
    startEndMat=cat(2,year(startNStime)',month(startNStime)',day(startNStime)',...
        hour(startNStime)',minute(startNStime)',floor(second(startNStime)')-1,...
        year(endNStime)',month(endNStime)',day(endNStime)',...
        hour(endNStime)',minute(endNStime)',ceil(second(endNStime)')+1,repmat(1,length(startNStime),1));
    
    disp([num2str(size(startEndMat,1)),' nsCal events found.']);
    
    startEndAll=cat(1,startEndAll,startEndMat);
end