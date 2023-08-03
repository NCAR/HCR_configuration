% Write topo data
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc0, qc1, qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5'; % ecmwf or era5

formatOut = 'yyyymmdd_HHMM';

[modeldir outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

topodir=topoDir(project);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

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
    
    data.dummy=[];
    
    dataVars=fieldnames(data);
    
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    %% We have some data where lat/lon are zero
    lonZero=find(data.longitude==0);
    if ~isempty(lonZero)
        warning([num2str(length(lonZero)),' zero longitudes replaced.']);
    end
    data.longitude(lonZero)=data.longitude(lonZero-1);

    latZero=find(data.latitude==0);
    if ~isempty(latZero)
        warning([num2str(length(latZero)),' zero latitudes replaced.']);
    end
    data.latitude(latZero)=data.latitude(latZero-1);
    
    %% Topo data
    disp('Getting topo data ...');
    modelData=[];
    [modelData.topo modelData.topolon modelData.topolat]=read_gtopo30(topodir,data.longitude,data.latitude);
    
    %% Interpolate
    disp('Interpolating ...');
    
    % Topo
    surfData.zHCR=interpn(modelData.topolon',modelData.topolat',modelData.topo',...
        wrapTo360(data.longitude),data.latitude);
    
    %% Save
    disp(['Saving topo data ...']);
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'timeHCR');
    topo=surfData.zHCR;
    save([outdir,whichModel,'.topo.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(ii),'.mat'],'topo');    
end

