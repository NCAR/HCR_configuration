% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.2';

plotInds=0;
% plotInds=(1:50:500);

outTime=0.1; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.1; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'airMotion/cases/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airMotion_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aircraft files
flightDir=[dataDirTS(1:end-27),'GV/highRate/'];
if ~exist(flightDir, 'dir')
    flightDir=[dataDirTS(1:end-27),'GV/lowRate/'];
end
flightFilesAll=dir([flightDir,'*.nc']);
flightsFile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];
flightsList=table2array(readtable(flightsFile));
flightStarts=datetime(flightsList(:,1:6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

% For FIR filter
Fnorm=10/(2000); % Normalized frequency
firFilt=designfilt("lowpassfir",FilterOrder=70,CutoffFrequency=Fnorm);
filtShift=mean(grpdelay(firFilt));

for aa=1:length(caseStart)
    tic

    plotTimeAll=[];
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% GV data
    flightDiff=flightStarts-startTime;
    gvInd=find(flightDiff>0,1)-1;
    gvFile=[flightDir,flightFilesAll(gvInd).name];
    aircraftVel=ncread(gvFile,'WIC');
    aircraftAlt=ncread(gvFile,'GGALT');
    aircraftTimeIn=ncread(gvFile,'Time');
    aircraftTime=datetime(flightStarts(gvInd).Year,flightStarts(gvInd).Month,flightStarts(gvInd).Day)+seconds(aircraftTimeIn);
    if min(size(aircraftVel))~=1
        aircraftTimeHR=repmat(aircraftTime',25,1);
        addSecs=0:1/25:1;
        addSecs(end)=[];
        addSecs=repmat(addSecs',1,length(aircraftTime));
        aircraftTimeHR=aircraftTimeHR+seconds(addSecs);

        aircraftDataHR=timetable(aircraftTimeHR(:),aircraftVel(:));
        aircraftDataLR=timetable(aircraftTime,aircraftAlt);

        aircraftData=synchronize(aircraftDataHR,aircraftDataLR,'first','linear');
    else
        if min(size(aircraftAlt))~=1
            aircraftAlt=aircraftAlt(1,:);
        end
        aircraftData=timetable(aircraftTime(:),aircraftVel(:),aircraftAlt(:));
    end

    %% CfRadial Moments
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.DBZ=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_CORR=[];
    dataCF.VEL_MASKED=[];
    dataCF.MELTING_LAYER=[];
    dataCF.ECHO_TYPE_2D=[];
    dataCF.CONVECTIVITY=[];
    dataCF.eastward_velocity=[];
    dataCF.northward_velocity=[];
         
    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);
    dataCF.beamWidth=ncread(fileListMoments{1},'radar_beam_width_v');

    % Find width correction
    velAircraft=sqrt(dataCF.eastward_velocity.^2+dataCF.northward_velocity.^2);
    widthCorrDelta=abs(0.3.*velAircraft.*sin(deg2rad(dataCF.elevation)).*deg2rad(dataCF.beamWidth));
    
    % Velocity bias term
    velBiasCorrection=dataCF.VEL_CORR-dataCF.VEL_RAW;
        
    %% Time series
    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    if isempty(fileListTS)
        warning('No data found.');
        continue
    end

    for bb=1:length(fileListTS)

        if bb==1
            data=[];
            data.IVc=[];
            data.QVc=[];
           
            vars=fieldnames(data);

            %% Load data TS
            disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);
           
            data=read_TsArchive_iwrf_bulk(fileListTS{1},data);

            % Find available times
            timeTest=data.time';
            timeTest(data.time<startTime | data.time>endTime)=[];
                                  
            % Set up de-aliasing
            toVel=data.lambda/(mode(data.prt)*4);

            % De-alias bias correction
            checkFold=[2,4,6];

            deAliasMaskB=zeros(size(velBiasCorrection));
            for jj=1:3
                deAliasMaskB(velBiasCorrection>checkFold(jj)*toVel-5)=checkFold(jj)*toVel;
                deAliasMaskB(velBiasCorrection<-(checkFold(jj)*toVel-5))=-checkFold(jj)*toVel;
            end
            velBiasCorrection=velBiasCorrection-deAliasMaskB;
        else
            %% Trimm first file
            data=trimFirstFile(data,endInd);

            goodTimes(1:end-1)=[];

            %% Load data TS
            disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);

            data=read_TsArchive_iwrf_bulk(fileListTS{bb},data);

            % Find available times
            timeTest=data.time';
            timeTest(data.time<goodTimes | data.time>endTime)=[];
        end

        TTdata=timetable(timeTest,ones(size(timeTest)));
        synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(outTime));
        goodTimes=synchTT.timeTest(synchTT.Var1>0);

        if bb<length(fileListTS)
            beamNum=length(goodTimes)-1;
        else
            beamNum=length(goodTimes);
        end

        %% Calculate moments

        disp('Calculating moments ...')
        momentsTimeOne.powerV=nan(size(data.range,1),beamNum);
        momentsTimeOne.powerH=nan(size(data.range,1),beamNum);
        momentsTimeOne.velRaw=nan(size(data.range,1),beamNum);
        momentsTimeOne.velRawDeAliased=nan(size(data.range,1),beamNum);
        momentsTimeOne.vel=nan(size(data.range,1),beamNum);
        momentsTimeOne.width=nan(size(data.range,1),beamNum);
        momentsTimeOne.dbz=nan(size(data.range,1),beamNum);
        momentsTimeOne.snr=nan(size(data.range,1),beamNum);
        momentsTimeOne.skew=nan(size(data.range,1),beamNum);
        momentsTimeOne.kurt=nan(size(data.range,1),beamNum);
        momentsTimeOne.ldr=nan(size(data.range,1),beamNum);
        momentsTimeOne.ncp=nan(size(data.range,1),beamNum);
        momentsTimeOne.range=nan(size(data.range,1),beamNum);
        momentsTimeOne.altitude=nan(1,beamNum);
        momentsTimeOne.elevation=nan(1,beamNum);
        momentsTimeOne.azimuth_vc=nan(1,beamNum);
        momentsTimeOne.time=goodTimes(1:beamNum)';

        momentsSpecOne=momentsTimeOne;

        lowShoulderVelOne=single(nan(size(data.range,1),beamNum,1));
        highShoulderVelOne=single(nan(size(data.range,1),beamNum,1));

        peakVelsOne=single(nan(size(data.range,1),beamNum,2));
        peakPowsOne=single(nan(size(data.range,1),beamNum,2));

        % Loop through beams
        for ii=1:beamNum
            
            % Find start and end indices for beam
            [~,startInd]=min(abs(goodTimes(ii)-seconds(sampleTime/2)-data.time));
            [~,endInd]=min(abs(goodTimes(ii)+seconds(sampleTime/2)-data.time));

            sampleNum=endInd-startInd+1;

            % Window
            win=window(@hamming,sampleNum);  % Default window is Hamming
            winWeight=sampleNum/sum(win);
            winNorm=win*winWeight;

            % Trim data down to current beam
            dataThis=trimData(data,startInd,endInd);

            % IQ
            cIQ.v=winNorm'.*(dataThis.IVc+i*dataThis.QVc)./sqrt(sampleNum);
                        
            %% Other variables
            momentsTimeOne.range(:,ii)=dataThis.range;
            momentsTimeOne.elevation(ii)=median(dataThis.elevation);
            momentsTimeOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsTimeOne.altitude(ii)=median(dataThis.altitude);

            momentsSpecOne.range(:,ii)=dataThis.range;
            momentsSpecOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecOne.altitude(:,ii)=median(dataThis.altitude);

            %% Find time index in CF moments
            cfInd=find(abs(etime(datevec(dataCF.time),datevec(momentsTimeOne.time(ii))))<0.0001);

            if isempty(cfInd)
                continue
            end

            %% Time moments
            momentsTimeOne=calcMomentsTime(cIQ,ii,momentsTimeOne,dataThis);

            % Censor on CF moments
            momentsTimeOne.velRaw(isnan(dataCF.VEL_MASKED(:,cfInd)),ii)=nan;

            %% Correct time domain velocity folding

            deAliasDiffT=momentsTimeOne.velRaw(:,ii)+velBiasCorrection(:,cfInd)-dataCF.VEL_MASKED(:,cfInd);

            deAliasMaskT=zeros(size(deAliasDiffT));
           
            for jj=1:3
                deAliasMaskT(deAliasDiffT>checkFold(jj)*toVel-5)=checkFold(jj)*toVel;
                deAliasMaskT(deAliasDiffT<-(checkFold(jj)*toVel-5))=-checkFold(jj)*toVel;
            end

            momentsTimeOne.velRawDeAliased(:,ii)=momentsTimeOne.velRaw(:,ii)-deAliasMaskT;
            
            %% Correct time domain velocity for aircraft motion and bias
            momentsTimeOne.vel(:,ii)=momentsTimeOne.velRawDeAliased(:,ii)+velBiasCorrection(:,cfInd);

            %% Spectra
            [specPowerLin,specPowerDB]=getSpectra(cIQ);
          
             % Censor
            specPowerLin.V(isnan(momentsTimeOne.vel(:,ii)),:)=nan;
            specPowerDB.V(isnan(momentsTimeOne.vel(:,ii)),:)=nan;

            %% Remove noise, find edge points and peaks
            if ismember(ii,plotInds)
                plotTime=momentsTimeOne.time(ii);
                plotTimeAll=cat(1,plotTimeAll,plotTime);
            else
                plotTime=[];
            end

            % This step removes the noise, de-aliases, (and corrects for
            % spectral broadening)
            [powerRMnoiseDBsmooth,specVelAdj,peakVels,peakPows]=noisePeaksAirVel(specPowerDB.V, ...
                momentsTimeOne.velRawDeAliased(:,ii),dataThis,firFilt,filtShift,widthCorrDelta(cfInd),plotTime);
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

            % Remove aircraft motion
            specVelAdj=specVelAdj+velBiasCorrection(:,cfInd);
            specVelRMnoise=specVelRMnoise+velBiasCorrection(:,cfInd);
            peakVels=peakVels+velBiasCorrection(:,cfInd);

            %% Spectral moments

            momentsSpecOne=calcMomentsSpec_higherMoments(powerRMnoiseDBsmooth,specVelRMnoise,ii,momentsSpecOne,dataThis);
          
            %% Velocities at peaks and shoulders
            peakVelsOne(:,ii,:)=peakVels;
            peakPowsOne(:,ii,:)=peakPows;
           
            % Find shoulders
            lowShoulderVelOne(:,ii)=specVelRMnoise(:,1);            

            flipSpec=fliplr(specVelRMnoise);
            [~,highInds]=max(~isnan(flipSpec),[],2);
            highIndsLin=sub2ind(size(specVelRMnoise),1:size(specVelRMnoise,1),highInds');
            highShoulderVelOne(:,ii)=flipSpec(highIndsLin);

            % % Width correction of shoulders
            % lowShoulderVelOne(:,ii)=lowShoulderVelOne(:,ii)+widthCorrDelta(cfInd)/2;
            % highShoulderVelOne(:,ii)=highShoulderVelOne(:,ii)-widthCorrDelta(cfInd)/2;

        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            momentsSpec=momentsSpecOne;
          
            lowShoulderVel=lowShoulderVelOne;
            highShoulderVel=highShoulderVelOne;
            peakVelsAll=peakVelsOne;
            peakPowsAll=peakPowsOne;

        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            dataFields1=fields(momentsSpec);

            for hh=1:length(dataFields1)
                momentsSpec.(dataFields1{hh})=cat(2,momentsSpec.(dataFields1{hh}),momentsSpecOne.(dataFields1{hh}));
            end

            lowShoulderVel=cat(2,lowShoulderVel,lowShoulderVelOne);
            highShoulderVel=cat(2,highShoulderVel,highShoulderVelOne);
            peakVelsAll=cat(2,peakVelsAll,peakVelsOne);
            peakPowsAll=cat(2,peakPowsAll,peakPowsOne);
        end
       
    end
    
    %% Reverse up pointing direction
    shoulderLowVel=nan(size(lowShoulderVel));
    shoulderLowVel(:,momentsTime.elevation<=0)=lowShoulderVel(:,momentsTime.elevation<=0);
    shoulderLowVel(:,momentsTime.elevation>0)=-highShoulderVel(:,momentsTime.elevation>0);
    shoulderHighVel=nan(size(highShoulderVel));
    shoulderHighVel(:,momentsTime.elevation<=0)=highShoulderVel(:,momentsTime.elevation<=0);
    shoulderHighVel(:,momentsTime.elevation>0)=-lowShoulderVel(:,momentsTime.elevation>0);

    dataCF.VEL_MASKED(:,dataCF.elevation>0)=-dataCF.VEL_MASKED(:,dataCF.elevation>0);
    momentsSpec.velRaw(:,momentsSpec.elevation>0)=-momentsSpec.velRaw(:,momentsSpec.elevation>0);
    momentsSpec.skew(:,momentsSpec.elevation>0)=-momentsSpec.skew(:,momentsSpec.elevation>0);
    
    peakLowVel=nan(size(lowShoulderVel));
    peakLowVel(:,momentsTime.elevation<=0)=peakVelsAll(:,momentsTime.elevation<=0,1);
    peakLowVel(:,momentsTime.elevation>0)=-peakVelsAll(:,momentsTime.elevation>0,2);
    peakHighVel=nan(size(highShoulderVel));
    peakHighVel(:,momentsTime.elevation<=0)=peakVelsAll(:,momentsTime.elevation<=0,2);
    peakHighVel(:,momentsTime.elevation>0)=-peakVelsAll(:,momentsTime.elevation>0,1);

    peakLowPow=nan(size(lowShoulderVel));
    peakLowPow(:,momentsTime.elevation<=0)=peakPowsAll(:,momentsTime.elevation<=0,1);
    peakLowPow(:,momentsTime.elevation>0)=peakPowsAll(:,momentsTime.elevation>0,2);
    peakHighPow=nan(size(highShoulderVel));
    peakHighPow(:,momentsTime.elevation<=0)=peakPowsAll(:,momentsTime.elevation<=0,2);
    peakHighPow(:,momentsTime.elevation>0)=peakPowsAll(:,momentsTime.elevation>0,1);

    %% Clean up peak layers

    % Velocity
    peakLowVel(isnan(peakHighVel))=nan;
    peakHighVel(isnan(peakLowVel))=nan;

    densMask=nan(size(peakLowVel));
    densMask(~isnan(peakLowVel))=1;

    [~,multiDens]=fast_nd_mean(densMask,[11,11]);

    multiMask=multiDens>15;
    multiMask=imclose(multiMask,strel('disk',3));
    multiMask=imopen(multiMask,strel('disk',1));
    multiMask=bwareaopen(multiMask,51);

    peakLowVel(multiMask==0)=nan;
    peakHighVel(multiMask==0)=nan;

    % Power
    peakLowPow(isnan(peakHighPow))=nan;
    peakHighPow(isnan(peakLowPow))=nan;

    peakLowPow(multiMask==0)=nan;
    peakHighPow(multiMask==0)=nan;

    %% Take time

    eSecs=toc;

    eData=momentsTime.time(end)-momentsTime.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Aircraft vel
     TTmoments=timetable(momentsTime.time',ones(size(momentsTime.time')));
     TTaircraft=synchronize(TTmoments,aircraftData,'first','nearest');
     TTaircraft.Properties.VariableNames={'Dummy';'Vel';'Alt'};
     TTaircraft.Properties.DimensionNames{1}='Time';

     %% Pull out data for comparison
     TTcompAC=timetable(momentsSpec.time',momentsSpec.velRaw(18,:)',shoulderLowVel(18,:)',momentsTime.dbz(18,:)',TTaircraft.Vel,TTaircraft.Alt, ...
         'VariableNames',{'totVel','lowVel','refl','aircraftVel','aircraftAlt'});
            
     cloudInds=find(TTcompAC.refl>-27 & TTcompAC.refl<0);
     precipInds=find(TTcompAC.refl>10);
     TTcloud=TTcompAC(cloudInds,:);
     TTprecip=TTcompAC(precipInds,:);

     TTcompAC.velCombined=nan(size(TTcompAC.Time));
     TTcompAC.velCombined(cloudInds)=TTcompAC.totVel(cloudInds);
     TTcompAC.velCombined(precipInds)=TTcompAC.lowVel(precipInds);

     cloudInds2D=find(momentsTime.dbz>-27 & momentsTime.dbz<0);
     precipInds2D=find(momentsTime.dbz>10);

     airMotion=nan(size(momentsTime.dbz));
     airMotion(cloudInds2D)=momentsSpec.velRaw(cloudInds2D);
     airMotion(precipInds2D)=shoulderLowVel(precipInds2D);

    %% Plot

    close all

    disp('Plotting velocities ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);
    momentsSpec.asl=HCRrange2asl(momentsSpec.range,momentsSpec.elevation,momentsSpec.altitude);

    plotAirMotion(momentsTime,momentsSpec,dataCF,shoulderLowVel,shoulderHighVel,peakLowVel,peakHighVel,TTaircraft,figdir,project,showPlot,plotTimeAll);
    
    %% Plot air motion
    close all
    plotAirMotionResults(momentsTime,airMotion,TTcompAC,TTcloud,TTprecip,figdir,project,showPlot);
end
