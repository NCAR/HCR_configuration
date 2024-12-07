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

figdir=[dataDirTS,'figsTS_rmNoise/'];

showPlot='on';
ylimUpper=7;
ylimLower=-0.1;
plotSpectra=1;

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/calcMoments_',project,'.txt'];
if plotSpectra
    sTimeFile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/spectraTimes_',project,'.txt'];
end

freqStr=strfind(freqData,'hz');
outFreq=str2num(freqData(1:freqStr-1)); % Desired output frequency in Hz
timeSpan=1/outFreq;

plotTimes=[];
plotRangeKM=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotSpectra
    plotTimeRange=readtable(sTimeFile);
    plotTimes=datetime(plotTimeRange.Var1,plotTimeRange.Var2,plotTimeRange.Var3, ...
        plotTimeRange.Var4,plotTimeRange.Var5,plotTimeRange.Var6);
    plotRangeKM=plotTimeRange.Var7;
end

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

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
    data.eastward_velocity=[];
    data.northward_velocity=[];
    data.vertical_velocity=[];
    data.azimuth_vc=[];

    data=readHCRts(fileListTS,data,startTime-seconds(timeSpan),endTime+seconds(timeSpan));

    %% Calculate moments

    disp('Calculating moments ...')

    % Find available times
    timeRound=dateshift(data.time,'start','minute')+seconds(round(second(data.time),1));
    timeRound=unique(timeRound);
    goodTimes=timeRound(timeRound>=startTime & timeRound<=endTime);

    beamNum=length(goodTimes);

    momentsTime.powerV=nan(size(data.range,1),beamNum);
    momentsTime.powerH=nan(size(data.range,1),beamNum);
    momentsTime.velRaw=nan(size(data.range,1),beamNum);
    momentsTime.vel=nan(size(data.range,1),beamNum);
    momentsTime.width=nan(size(data.range,1),beamNum);
    momentsTime.dbz=nan(size(data.range,1),beamNum);
    momentsTime.snr=nan(size(data.range,1),beamNum);
    momentsTime.skew=nan(size(data.range,1),beamNum);
    momentsTime.kurt=nan(size(data.range,1),beamNum);
    momentsTime.ldr=nan(size(data.range,1),beamNum);
    momentsTime.range=nan(size(data.range,1),beamNum);
    momentsTime.asl=nan(size(data.range,1),beamNum);
    momentsTime.elevation=nan(1,beamNum);
    momentsTime.eastward_velocity=nan(1,beamNum);
    momentsTime.northward_velocity=nan(1,beamNum);
    momentsTime.vertical_velocity=nan(1,beamNum);
    momentsTime.azimuth_vc=nan(1,beamNum);
    momentsTime.time=goodTimes;
       
    momentsSpec=momentsTime;
    momentsSpecRMnoiseSmooth=momentsTime;
    momentsSpecRMnoise=momentsTime;

  % Loop through beams
    for ii=1:beamNum

        % Find start and end indices for beam
        [~,startInd]=min(abs(etime(datevec(goodTimes(ii)-seconds(timeSpan/2)),datevec(data.time))));
        [~,endInd]=min(abs(etime(datevec(goodTimes(ii)+seconds(timeSpan/2)),datevec(data.time))));
        
        sampleNum=endInd-startInd+1;

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        % Trim data down to current beam
        dataThis=trimData(data,startInd,endInd);
        
        % IQ
        cIQ.v=winNorm'.*(dataThis.IVc+i*dataThis.QVc)./sqrt(sampleNum);
        cIQ.h=winNorm'.*(dataThis.IHx+i*dataThis.QHx)./sqrt(sampleNum);

        %% Other variables
        momentsTime.range(:,ii)=dataThis.range;
        momentsTime.elevation(ii)=median(dataThis.elevation);
        momentsTime.eastward_velocity(ii)=median(dataThis.eastward_velocity);
        momentsTime.northward_velocity(ii)=median(dataThis.northward_velocity);
        momentsTime.vertical_velocity(ii)=median(dataThis.vertical_velocity);
        momentsTime.azimuth_vc(ii)=median(dataThis.azimuth_vc);
        momentsTime.asl(:,ii)=median(dataThis.asl,2);

        momentsSpec.range(:,ii)=dataThis.range;
        momentsSpec.elevation(ii)=median(dataThis.elevation);
        momentsSpec.eastward_velocity(ii)=median(dataThis.eastward_velocity);
        momentsSpec.northward_velocity(ii)=median(dataThis.northward_velocity);
        momentsSpec.vertical_velocity(ii)=median(dataThis.vertical_velocity);
        momentsSpec.azimuth_vc(ii)=median(dataThis.azimuth_vc);
        momentsSpec.asl(:,ii)=median(dataThis.asl,2);

        momentsSpecRMnoiseSmooth.range(:,ii)=dataThis.range;
        momentsSpecRMnoiseSmooth.elevation(ii)=median(dataThis.elevation);
        momentsSpecRMnoiseSmooth.eastward_velocity(ii)=median(dataThis.eastward_velocity);
        momentsSpecRMnoiseSmooth.northward_velocity(ii)=median(dataThis.northward_velocity);
        momentsSpecRMnoiseSmooth.vertical_velocity(ii)=median(dataThis.vertical_velocity);
        momentsSpecRMnoiseSmooth.azimuth_vc(ii)=median(dataThis.azimuth_vc);
        momentsSpecRMnoiseSmooth.asl(:,ii)=median(dataThis.asl,2);

        momentsSpecRMnoise.range(:,ii)=dataThis.range;
        momentsSpecRMnoise.elevation(ii)=median(dataThis.elevation);
        momentsSpecRMnoise.eastward_velocity(ii)=median(dataThis.eastward_velocity);
        momentsSpecRMnoise.northward_velocity(ii)=median(dataThis.northward_velocity);
        momentsSpecRMnoise.vertical_velocity(ii)=median(dataThis.vertical_velocity);
        momentsSpecRMnoise.azimuth_vc(ii)=median(dataThis.azimuth_vc);
        momentsSpecRMnoise.asl(:,ii)=median(dataThis.asl,2);

        %% Time moments
        momentsTime=calcMomentsTime(cIQ,ii,momentsTime,dataThis);

        %% Spectral moments
        [specPowerLin,specPowerDB]=getSpectra(cIQ);

        % Move peak of spectra to middle
        [specPowerDBadj,specVelAdj]=adjSpecBoundsV(specPowerDB.V,momentsTime.velRaw(:,ii),dataThis);

        % Power fields
        momentsSpec=calcMomentsSpec_powerFields(specPowerLin,ii,momentsSpec,dataThis);

        % Higher order moments
        momentsSpec=calcMomentsSpec_higherMoments(specPowerDBadj,specVelAdj,ii,momentsSpec,dataThis);

        %% Remove noise
        [powerDBsmooth,powerRMnoiseDBsmooth]=rmNoiseSpec_10Hz(specPowerDBadj);
        powerRMnoiseDB=specPowerDBadj;
        powerRMnoiseDB(isnan(powerRMnoiseDBsmooth))=nan;
        specVelRMnoise=specVelAdj;
        specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

        % Power fields
        momentsSpecRMnoiseSmooth=calcMomentsSpec_powerFields(10.^(powerRMnoiseDBsmooth./10),ii,momentsSpecRMnoiseSmooth,dataThis);
        momentsSpecRMnoise=calcMomentsSpec_powerFields(10.^(powerRMnoiseDB./10),ii,momentsSpecRMnoise,dataThis);

        % Higher order moments
        momentsSpecRMnoiseSmooth=calcMomentsSpec_higherMoments(powerRMnoiseDBsmooth,specVelRMnoise,ii,momentsSpecRMnoiseSmooth,dataThis);
        momentsSpecRMnoise=calcMomentsSpec_higherMoments(powerRMnoiseDB,specVelRMnoise,ii,momentsSpecRMnoise,dataThis);

        %% Plot spectra
        if plotSpectra
            plotTimeDiff=abs(etime(datevec(goodTimes(ii)),datevec(plotTimes)));
            plotInd=find(plotTimeDiff<0.04);
            if ~isempty(plotInd)
                disp('Plotting spectra ...')
                close all
                plotSpectraExamplesRMnoise(data,momentsTime,momentsSpec,momentsSpecRMnoise,momentsSpecRMnoiseSmooth, ...
                    specVelAdj,specPowerDBadj,powerDBsmooth,powerRMnoiseDB,powerRMnoiseDBsmooth, ...
                    plotRangeKM,plotInd,ii,ylimUpper,figdir,project);                
            end
        end

    end

    %% Plot

    close all

    disp('Plotting moments ...');

    plotMomentsCompare(momentsTime,figdir,project,'Time',ylimUpper,showPlot,plotTimes,plotRangeKM);
    plotMomentsCompare(momentsTime,figdir,project,'Time',ylimUpper,showPlot,plotTimes,plotRangeKM);
    plotMomentsCompare(momentsSpec,figdir,project,'Spec',ylimUpper,showPlot,plotTimes,plotRangeKM);
    plotMomentsCompare(momentsSpecRMnoise,figdir,project,'SpecRMnoise',ylimUpper,showPlot,plotTimes,plotRangeKM);
    plotMomentsCompare(momentsSpecRMnoiseSmooth,figdir,project,'SpecRMnoiseSmooth',ylimUpper,showPlot,plotTimes,plotRangeKM);

    plotMomentsDiff(momentsTime,momentsSpecRMnoise,figdir,project,ylimUpper,showPlot);
end
