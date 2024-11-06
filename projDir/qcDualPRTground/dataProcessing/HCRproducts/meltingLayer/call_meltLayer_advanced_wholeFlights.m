% Analyze HCR clouds

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';
whichModel='era5';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

saveData=1;

plotYes=1;
showPlot='off';

thresholds.meltProbLow=0.4;
thresholds.meltProbHigh=0.55;

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'meltLayer/wholeIOPs/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% Loop through cases

for aa=1:size(caseList,1)

    disp(['IOP ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    %% Load data

    disp('Loading data ...');

    data=[];

    data.DBZ=[];
    data.VEL=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG_long=[];
    
    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
    end

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.FLAG=data.FLAG_long;
    data=rmfield(data,'FLAG_long');

    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end

    % % SPICULE has noisy LDR data that needs to be pre-processed
    % if strcmp(project,'spicule')
    %     disp('Pre-processing LDR ...');
    %     data.LDR=preProcessLDR(data.LDR);
    % end

    %% Find melting layer

    disp('Finding melting layer ...')
    data=f_meltLayer_advanced(data,thresholds,figdir);

    %% Save
    if saveData
        disp('Saving meltLayer field ...')

        meltLayer=data.meltLayer;
        save([outdir,whichModel,'.meltLayer.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'meltLayer');
    end

    %% Plot in hourly increments

    if plotYes
        disp('Plotting ...');
        ylimits=[0,6];
        startPlot=startTime;

        while startPlot<endTime

            endPlot=startPlot+minutes(20);
            indsTest=find(data.time>=startPlot & data.time<=endPlot);
            if length(indsTest)==0
                startPlot=endPlot;
                continue
            end
            
            newInds=indsTest(1):round(length(indsTest)/2000):indsTest(end);

            % Resample for plotting
            newDBZ=data.DBZ(:,newInds);
            if sum(sum(~isnan(newDBZ)))>300

                newLDR=data.LDR(:,newInds);
                velPlot=data.VEL;
                velPlot(:,data.elevation>0)=-velPlot(:,data.elevation>0);
                newVEL=velPlot(:,newInds);
                newVELdiff=data.velDiff(:,newInds);
                newDBZdiff=data.dbzDiff(:,newInds);
                newMeltLayer=data.meltLayer(:,newInds);
                newASL=data.asl(:,newInds);
                newTime=data.time(newInds);

                newProb=data.meltProb(:,newInds);
                newProb(newProb<0.1)=nan;

                timeInds=find(data.time>=newTime(1) & data.time<=newTime(end));
                timeForMask=data.time(timeInds);
                aslForMask=data.asl(:,timeInds);
                maskForPlot=data.meltMask(:,timeInds);
                timeMat=repmat(timeForMask,size(data.DBZ,1),1);

                close all

                meltTestPlot1;
                meltTestPlot2;
            end
            startPlot=endPlot;
        end
    end
end