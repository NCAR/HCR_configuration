% Call cloud classification script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';
whichModel='era5';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

saveTime=0;

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    disp('Loading data ...');
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.ECHO_TYPE_2D=[];
    data.FLAG=[];
    data.TEMP=[];
    data.TOPO=[];
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
 
    ylimUpper=(max(data.asl(~isnan(data.ECHO_TYPE_2D)))./1000)+0.5;
   
    %% Join clouds that are separated by missing

    missingData=any(data.FLAG>=10,1);
    joinedEcho=data.ECHO_TYPE_2D(:,missingData==0);
    joinedTemp=data.TEMP(:,missingData==0);
    joinedAsl=data.asl(:,missingData==0);
    joinedElev=data.elevation(missingData==0);
    joinedTopo=data.TOPO(missingData==0);

    %% Create cloudID

    disp('Creating cloud ID ...')

    minCloudSizePix=1000;

    cloudID=makeCloudID(joinedEcho,minCloudSizePix);

    %% Cloud classification

    disp('Finding cloud class ...')

    cloudClassJoined=findCloudClass(joinedEcho,cloudID,joinedTemp,joinedElev,joinedTopo,joinedAsl);

    cloudClass=nan(size(data.ECHO_TYPE_2D));
    cloudClass(:,missingData==0)=cloudClassJoined;

    % Fill in small with not classified
    cloudClass(~isnan(data.ECHO_TYPE_2D) & isnan(cloudClass))=0;

     %% Save
    disp('Saving cloudClass field ...')
    
    save([outdir,whichModel,'.cloudClass.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'cloudClass');
        
    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end
       
end