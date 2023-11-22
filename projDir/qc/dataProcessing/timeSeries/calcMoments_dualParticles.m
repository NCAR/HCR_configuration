% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.1';

dataDirTS=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDirTS,'figsTS_dualParticles/'];

showPlot='on';
ylimUpper=5.2;
ylimLower=-0.1;

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/dualParticles_',project,'.txt'];

freqStr=strfind(freqData,'hz');
outFreq=str2num(freqData(1:freqStr-1)); % Desired output frequency in Hz
timeSpan=1/outFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Load data TS
    disp("Getting time series data ...");

    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    data=[];
    data.IVc=[];
    data.QVc=[];
    data.IHx=[];
    data.QHx=[];

    data=readHCRts(fileListTS,data,startTime,endTime);

    %% Calculate moments

    disp('Calculating moments ...')

    beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

    momentsTime.powerV=nan(size(data.range,1),beamNum);
    momentsTime.powerH=nan(size(data.range,1),beamNum);
    momentsTime.vel=nan(size(data.range,1),beamNum);
    momentsTime.width=nan(size(data.range,1),beamNum);
    momentsTime.dbz=nan(size(data.range,1),beamNum);
    momentsTime.snr=nan(size(data.range,1),beamNum);
    momentsTime.skew=nan(size(data.range,1),beamNum);
    momentsTime.kurt=nan(size(data.range,1),beamNum);
    momentsTime.ldr=nan(size(data.range,1),beamNum);

    momentsSpecRMnoise=momentsTime;

    momentsVelDual=nan(size(data.range,1),beamNum,2);

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

        cIQ.v=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);
        cIQ.h=winNorm'.*(data.IHx(:,startInd:endInd)+i*data.QHx(:,startInd:endInd))./sqrt(sampleNum);

        data.prtThis=data.prt(startInd:endInd);

        %% Time moments
        momentsTime=calcMomentsTime(cIQ,ii,momentsTime,data);

        %% Spectral moments
        [specPowerLin,specPowerDB]=getSpectra(cIQ);
        
        % Move peak of spectra to middle
        [specPowerDBadj,specVelAdj]=adjSpecBoundsV(specPowerDB.V,momentsTime.vel(:,ii),sampleNum,data);

        %% Remove noise
        [powerDBsmooth,powerRMnoiseDBsmooth]=rmNoiseSpec(specPowerDBadj,specVelAdj,sampleNum);
        powerRMnoiseDB=specPowerDBadj;
        powerRMnoiseDB(isnan(powerRMnoiseDBsmooth))=nan;
        specVelRMnoise=specVelAdj;
        specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

        % Power fields
        momentsSpecRMnoise=calcMomentsSpec_powerFields(10.^(powerRMnoiseDB./10),ii,momentsSpecRMnoise,data);

        % Higher order moments
        momentsSpecRMnoise=calcMomentsSpec_higherMoments(powerRMnoiseDB,specVelRMnoise,ii,momentsSpecRMnoise,data);

        %% Find regions with dual particle species

        momentsVelDual=findDualParticles_test(powerRMnoiseDBsmooth,specVelRMnoise,specPowerDBadj,momentsVelDual,ii);

        %% Other processing
        
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

    %% Plot

    close all

    disp('Plotting velocities ...');

    plotVelocities(data,momentsVelDual,momentsTime,timeBeams,figdir,project,ylimUpper,flipYes,showPlot);

end
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');