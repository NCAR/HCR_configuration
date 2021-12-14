% Remove stripes in HX
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc0, qc1, qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, or 2hz
%whichModel='narr'; % ecmwf or era5

formatOut = 'yyyymmdd_HHMM';

%[modeldir outdir]=modelDir(project,whichModel,freqData);

topodir=topoDir(project);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

startEndTimesAll=[];

%% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    %% HCR data
    disp('Getting HCR data ...');
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    data.DBMHX=[];
    
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
    
    %% Find indices where correction needs to be applied
    medianHXend=median(data.DBMHX(760:770,:));
    
    corrInd=medianHXend>-95;
    
    if max(corrInd>0)
        disp('High DBMHX found.');
    end
    
    % Find time indices
    corrIndFill=nan(size(corrInd));
    corrIndFill(corrInd==1)=1;
    
    corrIndFill=movmean(corrIndFill,1800,'omitnan');
    corrIndFill(isnan(corrIndFill))=0;
    corrIndDiff=diff(corrIndFill);
    
    startTimes=data.time(corrIndDiff==1);
    endTimes=data.time(corrIndDiff==-1);

    if isempty(startTimes) & isempty(endTimes)
        continue
    end
    
    if startTimes(1)>endTimes(1)
        startTimes=cat(2,data.time(1),startTimes);
    end
    if length(startTimes)~=length(endTimes)
        endTimes=cat(2,endTimes,data.time(end));
    end
    
    startEndTimes=cat(2,startTimes',endTimes');
    startEndTimesAll=cat(1,startEndTimesAll,startEndTimes);
    
    disp(startEndTimes);
    
end
%% Save
outTable=table(startEndTimesAll);
outTable.startEndTimesAll.Format = 'yyyy-MM-dd HH:mm:ss';
writetable(outTable,[project,'HighDBMHXtimes.txt'],'delimiter',' ');