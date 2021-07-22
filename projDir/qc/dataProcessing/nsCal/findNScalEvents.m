% Find nsCal events
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='field'; % field, qc0, qc1, qc2
qcVersion='v0.1';
freq='10hz';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data_TFRF.txt'];

caseList = table2array(readtable(infile));

%indir=HCRdir(project,quality,qcVersion,freq);
indir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc0/cfradial/moments/10hz/';

%% Run processing

startEndAll=[];

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    data.DBMVC=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isempty(data.time)
        disp('No data found.');
        startTime=endTime;
        continue
    end
    
    % Check if all variables were found
    for kk=1:length(dataVars)
        if ~isfield(data,dataVars{kk})
            dataVars{kk}=[];
        end
    end
    
    %% Find nsCal
    firstGate=data.DBMVC(1,:);
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