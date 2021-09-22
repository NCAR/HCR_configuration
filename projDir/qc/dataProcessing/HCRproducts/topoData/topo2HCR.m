% Write topo data
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc0, qc1, qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5'; % ecmwf or era5

formatOut = 'yyyymmdd_HHMM';

[modeldir outdir]=modelDir(project,whichModel,freqData);

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
    
    data.asl=HCRrange2asl(data.range,data.elevation,data.altitude);
    
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

