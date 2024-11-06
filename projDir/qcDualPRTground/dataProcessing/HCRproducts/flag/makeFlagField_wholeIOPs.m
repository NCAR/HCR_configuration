% Analyze HCR clouds

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='meow'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz_combined'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';
whichModel='era5';
ls='_long'; % _short or _long

saveTime=1;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

directories.dataDir=HCRdir(project,quality,qcVersion,freqData);

[~,directories.modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

outdir=directories.modeldir;

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

%% Load data
for mm=1:size(caseList,1)
    disp(['IOP ',num2str(mm)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(mm,1:6));
    endTime=datetime(caseList(mm,7:12));
    
    data=[];
    
    data.(['DBZ',ls])=[];
    data.(['WIDTH',ls])=[];
    data.(['DBMVC',ls])=[];
               
    % Make list of files within the specified time frame
    fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
      
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    %% Mask data
    
    disp('Making flag field.')
    
    [maskData,antStat]=echoMask(data,ls);
    
    
    %% Save
    disp('Saving flag field.')
    
    flagField=maskData;
    save([outdir,whichModel,'.flagfield',ls,'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(mm),'.mat'],'flagField','-v7.3');
    
    antennaStat=antStat;
    save([outdir,whichModel,'.antstat',ls,'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(mm),'.mat'],'antennaStat');
        
    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time',ls,'.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(mm),'.mat'],'timeHCR');
    end
end