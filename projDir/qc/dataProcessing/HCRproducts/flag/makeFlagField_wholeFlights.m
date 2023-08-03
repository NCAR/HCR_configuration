% Analyze HCR clouds

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v3.1';
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

directories.dataDir=HCRdir(project,quality,qcVersion,freqData);

[~,directories.modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

outdir=directories.modeldir;

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%% Load data

for mm=1:size(caseList,1)
    disp(['Flight ',num2str(mm)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(mm,1:6));
    endTime=datetime(caseList(mm,7:12));
    
    data=[];
    
    data.DBZ=[];
    data.WIDTH=[];
    data.DBMVC=[];
    data.TOPO=[];
           
    % Make list of files within the specified time frame
    fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
      
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    %% Mask data
    
    disp('Making flag field.')
    
    [maskData,antStat]=echoMask(data);
    
    
    %% Save
    disp('Saving flag field.')
    
    flagField=maskData;
    save([outdir,whichModel,'.flagfield.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(mm),'.mat'],'flagField');
    
    antennaStat=antStat;
    save([outdir,whichModel,'.antstat.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(mm),'.mat'],'antennaStat');
        
end