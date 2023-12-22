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

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Load data TS
    disp("Getting time series data ...");

    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    data=[];
   
    [data,fileStartAll,fileEndAll]=readHCRts(fileListTS,data,startTime-seconds(timeSpan),endTime+seconds(timeSpan),0);

    %% Calculate moments

    disp('Calculating moments ...')

    % Find available times
    timeTest=data.time';
    timeTest(data.time<startTime | data.time>endTime)=[];
    TTdata=timetable(timeTest,ones(size(timeTest)));
    synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(timeSpan));
    goodTimes=synchTT.timeTest(synchTT.Var1>0);

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
    momentsTime.altitude=nan(1,beamNum);
    momentsTime.elevation=nan(1,beamNum);
    momentsTime.eastward_velocity=nan(1,beamNum);
    momentsTime.northward_velocity=nan(1,beamNum);
    momentsTime.vertical_velocity=nan(1,beamNum);
    momentsTime.azimuth_vc=nan(1,beamNum);
    momentsTime.time=goodTimes;

    momentsVelDualRaw=single(nan(size(data.range,1),beamNum,1));

    tic
    % Loop through beams
    for ii=1:beamNum

        %disp(datestr(goodTimes(ii),'yyyymmdd_HHMMSS.FFF'));
        if ii==111
            stop1=1;
        end

        % Find start and end indices for beam
        [~,startInd]=min(abs(goodTimes(ii)-seconds(timeSpan/20)-data.time));
        [~,endInd]=min(abs(goodTimes(ii)+seconds(timeSpan/20)-data.time));
       
        sampleNum=endInd-startInd+1;

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        % Find correct files
        startTest=startInd-fileStartAll;
        startTest(startTest<0)=[];
        startFile=length(startTest);

        endTest=endInd-fileEndAll;
        endTest(endTest<0)=[];
        endFile=length(endTest)+1;

        fileListTSthis=fileListTS(:,startFile:endFile);
        readStart=startInd-fileStartAll(startFile);

        % Load data
        dataThis=[];
        dataThis.IVc=[];
        dataThis.QVc=[];
        dataThis.IHx=[];
        dataThis.QHx=[];
        dataThis.eastward_velocity=[];
        dataThis.northward_velocity=[];
        dataThis.vertical_velocity=[];
        dataThis.azimuth_vc=[];
        
        dataThis=readHCRts_startEnd(fileListTSthis,dataThis,readStart,sampleNum);
        
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
        momentsTime.altitude(ii)=median(dataThis.altitude);

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
        [powerDBsmooth,powerRMnoiseDBsmooth]=rmNoiseSpec(specPowerDBadj);
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
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Sort out vel dual
    momentsVelDual=sortDualParticles(momentsVelDualRaw,momentsTime);

    %% Plot

    close all

    disp('Plotting velocities ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);

    plotVelocities(momentsVelDual,momentsTime,figdir,project,showPlot);

end
