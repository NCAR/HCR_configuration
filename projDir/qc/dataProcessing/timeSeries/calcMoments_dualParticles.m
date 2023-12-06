% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz'; % !!!!!!!!! Must be equal or less than one second !!!!!!!!!!!!!
qcVersion='v1.2';

dataDirTS=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDirTS,'figsTS_dualParticles/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/dualParticles_',project,'.txt'];

freqStr=strfind(freqData,'hz');
outFreq=str2num(freqData(1:freqStr-1)); % Desired output frequency in Hz
timeSpan=1/outFreq; % Desired time resolution in seconds

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
    TTdata=timetable(data.time',ones(size(data.time))');
    synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(timeSpan));
    goodTimes=synchTT.Time(synchTT.Var1>0);

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
    momentsTime.ncp=nan(size(data.range,1),beamNum);
    momentsTime.range=nan(size(data.range,1),beamNum);
    momentsTime.asl=nan(size(data.range,1),beamNum);
    momentsTime.elevation=nan(1,beamNum);
    momentsTime.eastward_velocity=nan(1,beamNum);
    momentsTime.northward_velocity=nan(1,beamNum);
    momentsTime.vertical_velocity=nan(1,beamNum);
    momentsTime.azimuth_vc=nan(1,beamNum);
    momentsTime.time=goodTimes;

    momentsVelDualRaw=nan(size(data.range,1),beamNum,1);

    tic
    % Loop through beams
    for ii=1:beamNum

        %disp(datestr(goodTimes(ii),'yyyymmdd_HHMMSS.FFF'));

        % Find start and end indices for beam
        [~,startInd]=min(abs(goodTimes(ii)-seconds(timeSpan/20)-data.time));
        [~,endInd]=min(abs(goodTimes(ii)+seconds(timeSpan/20)-data.time));
        
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

        %% Time moments
        momentsTime=calcMomentsTime(cIQ,ii,momentsTime,dataThis);

        % Censor on SNR and NCP
        censorY=momentsTime.snr(:,ii)<0 & momentsTime.ncp(:,ii)<0.2;
        censorY(isnan(momentsTime.snr(:,ii)))=1;
        censorY=~censorY;
        censorY=double(censorY);
        censorY(censorY==0)=nan;
        censorY=movmedian(censorY,7,'includemissing');
        censorY=movmedian(censorY,7,'omitmissing');
        momentsTime.vel(isnan(censorY),ii)=nan;

        %% Spectral moments
        [specPowerLin,specPowerDB]=getSpectra(cIQ);
        
        % Move peak of spectra to middle
        [specPowerDBadj,specVelAdj]=adjSpecBoundsV(specPowerDB.V,momentsTime.velRaw(:,ii),dataThis);

        % Censor
        specPowerDBadj(isnan(momentsTime.vel(:,ii)),:)=nan;
        specVelAdj(isnan(momentsTime.vel(:,ii)),:)=nan;

        %% Remove noise
        [powerDBsmooth,powerRMnoiseDBsmooth]=rmNoiseSpec(specPowerDBadj,specVelAdj);
        powerRMnoiseDB=specPowerDBadj;
        powerRMnoiseDB(isnan(powerRMnoiseDBsmooth))=nan;
        specVelRMnoise=specVelAdj;
        specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

        %% Find regions with dual particle species

        momentsVelDualRaw=findDualParticles(powerRMnoiseDBsmooth,specVelRMnoise,specPowerDBadj,momentsVelDualRaw,ii);

    end
    eSecs=toc;

    eData=momentsTime.time(end)-momentsTime.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp([num2str(timePerMin),' minutes per data minute.']);

    %% Sort out vel dual
    momentsVelDual=sortDualParticles(momentsVelDualRaw,momentsTime);

    %% Plot

    close all

    disp('Plotting velocities ...');

    plotVelocities(momentsVelDual,momentsTime,figdir,project,showPlot);

end
