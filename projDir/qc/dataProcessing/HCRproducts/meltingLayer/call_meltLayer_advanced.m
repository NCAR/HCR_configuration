% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

plotYes=1;
showPlot='off';

thresholds.meltProbLow=0.4; %0.55
thresholds.meltProbHigh=0.55; %0.7

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/meltLayer_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'meltLayer/cases/'];

if ~exist(figdir,'dir')
    mkdir(figdir)
end

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
    data.DBZ_MASKED=[];
    data.VEL_MASKED=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];
        
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end

    data.LDR(data.FLAG~=1)=nan;

    %% Find melting layer

    disp('Finding melting layer ...')
    data=f_meltLayer_advanced(data,thresholds,figdir);

    %% Plot
    if plotYes
        disp('Plotting ...')
        ylimits=[0,5];

        newInds=1:round(length(data.time)/2000):length(data.time);

        % Resample for plotting
        newDBZ=data.DBZ_MASKED(:,newInds);
        newLDR=data.LDR(:,newInds);
        velPlot=data.VEL_MASKED;
        velPlot(:,data.elevation>0)=-velPlot(:,data.elevation>0);
        newVEL=velPlot(:,newInds);
        newVELdiff=data.velDiff(:,newInds);
        newDBZdiff=data.dbzDiff(:,newInds);
        newASL=data.asl(:,newInds);
        newTime=data.time(newInds);

        newProb=data.meltProb(:,newInds);
        newProb(newProb<0.1)=nan;

        timeInds=find(data.time>=newTime(1) & data.time<=newTime(end));
        timeForMask=data.time(timeInds);
        aslForMask=data.asl(:,timeInds);
        maskForPlot=data.meltMask(:,timeInds);
        timeMat=repmat(timeForMask,size(data.DBZ_MASKED,1),1);

        close all

        meltTestPlot1;
        meltTestPlot2;
    end
end