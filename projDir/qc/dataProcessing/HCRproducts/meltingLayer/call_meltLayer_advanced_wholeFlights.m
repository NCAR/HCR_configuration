% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

plotYes=1;
showPlot='off';

saveData=1;

thresholds.meltProbLow=0.4;
thresholds.meltProbHigh=0.55;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'meltLayer/wholeFlights/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% Loop through cases

for aa=1:size(caseList,1)

    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    %% Load data

    disp('Loading data ...');

    data=[];

    data.DBZ_MASKED=[];
    data.VEL_MASKED=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];

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

    % SPICULE has noisy LDR data that needs to be pre-processed
    if strcmp(project,'spicule')
        disp('Pre-processing LDR ...');
        data.LDR=preProcessLDR(data.LDR);
    end

    %% Find melting layer

    disp('Finding melting layer ...')
    data=f_meltLayer_advanced(data,thresholds,figdir);

    %% Save
    if saveData
        disp('Saving meltLayer field ...')

        meltLayer=data.meltLayer;
        save([outdir,whichModel,'.meltLayer.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'meltLayer');

        iceLev=data.iceLev;
        save([outdir,whichModel,'.iceLevel.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'iceLev');
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
            newDBZ=data.DBZ_MASKED(:,newInds);
            if sum(sum(~isnan(newDBZ)))>300

                newLDR=data.LDR(:,newInds);
                velPlot=data.VEL_MASKED;
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
                timeMat=repmat(timeForMask,size(data.DBZ_MASKED,1),1);

                close all

                meltTestPlot1;
                meltTestPlot2;
            end
            startPlot=endPlot;
        end
    end
end