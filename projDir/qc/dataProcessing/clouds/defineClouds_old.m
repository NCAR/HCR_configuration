% Analyze HCR clouds

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

figdir=['/h/eol/romatsch/hcrCalib/mask/',project,'/defineClouds/'];
metadir=['/scr/sci/romatsch/data/HCRprods/mask/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

directories.dataDir=HCRdir(project,quality,freqData);

%% Identify time intervals with clouds
load([metadir,'/socrates.cloudInds.20180115_224025_to_20180116_052609.Flight1.mat']);
load([metadir,'/socrates.timeFlight.20180115_224025_to_20180116_052609.Flight1.mat']);

diffInds=diff(cloudInds);

startAll=timeFlight(find(diffInds==1)+1);
endAll=timeFlight(find(diffInds==-1));

isInd=find(diffInds~=0);
if diffInds(isInd(1))==-1
    startAll=cat(2,timeFlight(1),startAll);
end
if diffInds(isInd(end))==1
    endAll=cat(2,endAll,timeFlight(end));
end

if length(startAll)~=length(endAll)
    disp('Something is wrong with the start and end times.')
    return
end

singleInds=find(startAll==endAll);
startAll(singleInds)=[];
endAll(singleInds)=[];
    
%% Go through flight hours
for ii=1:length(startAll)
    startTime=startAll(ii);
    endTime=endAll(ii);
    
    % Desired variables. The variable name comies after the . and must be spelled exactly
    % like in the CfRadial file
    if exist('data')
        clear data
    end
    
    %% Load data
    
    data.DBZ=[];
    %data.VEL=[];
    %data.VEL_RAW=[];
    %data.VEL_CORR=[];
    data.WIDTH=[];
    %data.WIDTH_CORR=[];
    %data.DBMVC=[];
    %data.SNR=[];
    %data.NCP=[];
    %data.LDR=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        startTime=endTime;
        continue
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
end