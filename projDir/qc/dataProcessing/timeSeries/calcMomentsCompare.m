% Analyze HCR time series
clear all;
close all;

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.1';

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirTS,'figsTS/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airVel_',project,'.txt'];

freqStr=strfind(freqData,'hz');
outFreq=str2num(freqData(1:freqStr-1)); % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=7.5;
ylimLower=-0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

for aa=4:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Load data TS
    disp("Getting time series data ...");

    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    data=[];
    data.IVc=[];
    data.QVc=[];

    data=readHCRts(fileListTS,data,startTime,endTime);

    %% Calculate moments
    beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

    momentsTime.powerV=nan(size(data.range,1),beamNum);
    momentsTime.vel=nan(size(data.range,1),beamNum);
    momentsTime.width=nan(size(data.range,1),beamNum);
    momentsTime.dbz=nan(size(data.range,1),beamNum);
    momentsTime.snr=nan(size(data.range,1),beamNum);
    momentsTime.skew=nan(size(data.range,1),beamNum);
    momentsTime.kurt=nan(size(data.range,1),beamNum);
    
    momentsSpec.powerV=nan(size(data.range,1),beamNum);
    momentsSpec.vel=nan(size(data.range,1),beamNum);
    momentsSpec.width=nan(size(data.range,1),beamNum);
    momentsSpec.dbz=nan(size(data.range,1),beamNum);
    momentsSpec.snr=nan(size(data.range,1),beamNum);
    momentsSpec.skew=nan(size(data.range,1),beamNum);
    momentsSpec.kurt=nan(size(data.range,1),beamNum);

    timeBeams=[];

    startInd=1;
    endInd=1;
    ii=1;

    % Loop through beams
    while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

        % Find start and end indices for beam
        timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
        [minDiff,endInd]=min(abs(timeDiff+timeSpan));

        sampleNum=endInd-startInd+1;

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

        prt=mode(data.prt);

        %% Spectral moments
        momentsSpec=calcMomentsSpec(cIQv,sampleNum,ii,momentsSpec,data);

        %% Time moments
        momentsTime=calcMomentsTime(cIQv,prt,ii,momentsTime,data);
        
        timeBeams=[timeBeams;data.time(startInd)];

        if data.elevation(startInd)<0
            flipYes=1;
        else
            flipYes=0;
        end

        %% Next beam
        startInd=endInd+1;
        ii=ii+1;
    end

    momentsSpec.vel=momentsSpec.vel.*data.lambda/(4*pi*prt);

    %% Plot

    close all

    disp('Plotting ...');

    plotMomentsCompare(data,momentsSpec,timeBeams,figdir,project,'Spec',ylimUpper,flipYes,showPlot);
    plotMomentsCompare(data,momentsTime,timeBeams,figdir,project,'Time',ylimUpper,flipYes,showPlot);
end