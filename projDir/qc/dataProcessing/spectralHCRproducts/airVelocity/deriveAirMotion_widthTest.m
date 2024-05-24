% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc3';
freqData='10hz';
qcVersion='v3.2';

plotInds=0;
%plotInds=(1:50:500);

outTime=0.1; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.1; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'airMotion/cases/width/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airMotion_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

for aa=1:length(caseStart)
    tic

    plotTimeAll=[];
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% CfRadial Moments
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.DBZ=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_CORR=[];
    dataCF.VEL_MASKED=[];
    dataCF.eastward_velocity=[];
    dataCF.northward_velocity=[];
         
    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);
    dataCF.beamWidth=ncread(fileListMoments{1},'radar_beam_width_v');

    % Find width correction factor
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
        momentsTimeOne.widthCorr=nan(size(data.range,1),beamNum);
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

        momentsSpecBasicOne=momentsTimeOne;
        momentsSpecNoNoiseOne=momentsTimeOne;
        momentsSpecSmoothOne=momentsTimeOne;
        momentsSpecCorrectedOne=momentsTimeOne;

        % Loop through beams
        for ii=1:beamNum
            
            % Find start and end indices for beam
            [~,startInd]=min(abs(goodTimes(ii)-seconds(sampleTime/2)-data.time));
            [~,endInd]=min(abs(goodTimes(ii)+seconds(sampleTime/2)-data.time));

            sampleNum=endInd-startInd+1;

            if sampleNum<30
                continue
            end

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

            momentsSpecBasicOne.range(:,ii)=dataThis.range;
            momentsSpecBasicOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecBasicOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecBasicOne.altitude(:,ii)=median(dataThis.altitude);

            momentsSpecNoNoiseOne.range(:,ii)=dataThis.range;
            momentsSpecNoNoiseOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecNoNoiseOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecNoNoiseOne.altitude(:,ii)=median(dataThis.altitude);

            momentsSpecSmoothOne.range(:,ii)=dataThis.range;
            momentsSpecSmoothOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecSmoothOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecSmoothOne.altitude(:,ii)=median(dataThis.altitude);

            momentsSpecCorrectedOne.range(:,ii)=dataThis.range;
            momentsSpecCorrectedOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecCorrectedOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecCorrectedOne.altitude(:,ii)=median(dataThis.altitude);

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

            %% Correct time domain width
            widthSquares=momentsTimeOne.width(:,ii).^2-widthCorrDelta(ii).^2;
            widthSquares(widthSquares<0.1)=0.01;
            momentsTimeOne.widthCorr(:,ii)=sqrt(widthSquares);

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
            [powerRMnoiseDBcorrected,powerRMnoise,powerRMnoiseDBsmooth,specVelAdj,specVelAdjSmooth]=noisePeaksAirVel_widthTest(specPowerDB.V, ...
                momentsTimeOne.velRawDeAliased(:,ii),dataThis,widthCorrDelta(cfInd));
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerRMnoiseDBcorrected))=nan;

            specVelRMnoiseSmooth=specVelAdjSmooth;
            specVelRMnoiseSmooth(isnan(powerRMnoiseDBsmooth))=nan;

            % Remove aircraft motion
            specVelAdj=specVelAdj+velBiasCorrection(:,cfInd);
            specVelRMnoise=specVelRMnoise+velBiasCorrection(:,cfInd);
            
            %% Spectral moments

            momentsSpecBasicOne=calcMomentsSpec_higherMoments(specPowerDB.V,specVelAdj,ii,momentsSpecBasicOne,dataThis);
            momentsSpecNoNoiseOne=calcMomentsSpec_higherMoments(powerRMnoise,specVelRMnoiseSmooth,ii,momentsSpecNoNoiseOne,dataThis);
            momentsSpecSmoothOne=calcMomentsSpec_higherMoments(powerRMnoiseDBsmooth,specVelRMnoiseSmooth,ii,momentsSpecSmoothOne,dataThis);
            momentsSpecCorrectedOne=calcMomentsSpec_higherMoments(powerRMnoiseDBcorrected,specVelRMnoise,ii,momentsSpecCorrectedOne,dataThis);
        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            momentsSpecBasic=momentsSpecBasicOne;
            momentsSpecNoNoise=momentsSpecNoNoiseOne;
            momentsSpecSmooth=momentsSpecSmoothOne;
            momentsSpecCorrected=momentsSpecCorrectedOne;
        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            dataFields1=fields(momentsSpecBasic);

            for hh=1:length(dataFields1)
                momentsSpecBasic.(dataFields1{hh})=cat(2,momentsSpecBasic.(dataFields1{hh}),momentsSpecBasicOne.(dataFields1{hh}));
                momentsSpecNoNoise.(dataFields1{hh})=cat(2,momentsSpecNoNoise.(dataFields1{hh}),momentsSpecNoNoiseOne.(dataFields1{hh}));
                momentsSpecSmooth.(dataFields1{hh})=cat(2,momentsSpecSmooth.(dataFields1{hh}),momentsSpecSmoothOne.(dataFields1{hh}));
                momentsSpecCorrected.(dataFields1{hh})=cat(2,momentsSpecCorrected.(dataFields1{hh}),momentsSpecCorrectedOne.(dataFields1{hh}));
            end
        end
       
    end
    
    %% Take time

    eSecs=toc;

    eData=momentsTime.time(end)-momentsTime.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Match times

    tsRound=dateshift(momentsTime.time,'start','minute')+seconds(round(second(momentsTime.time),1));
    cfRound=dateshift(dataCF.time,'start','minute')+seconds(round(second(dataCF.time),1));

    indTs=ismember(tsRound,cfRound);
    indCf=ismember(cfRound,tsRound);

    fieldsTs=fieldnames(momentsTime);
    for ll=1:length(fieldsTs)
        momentsTime.(fieldsTs{ll})=momentsTime.(fieldsTs{ll})(:,indTs);
        momentsSpecBasic.(fieldsTs{ll})=momentsSpecBasic.(fieldsTs{ll})(:,indTs);
        momentsSpecNoNoise.(fieldsTs{ll})=momentsSpecNoNoise.(fieldsTs{ll})(:,indTs);
        momentsSpecSmooth.(fieldsTs{ll})=momentsSpecSmooth.(fieldsTs{ll})(:,indTs);
        momentsSpecCorrected.(fieldsTs{ll})=momentsSpecCorrected.(fieldsTs{ll})(:,indTs);
    end

    dataCF=rmfield(dataCF,'beamWidth');

    fieldsCf=fieldnames(dataCF);
    for ll=1:length(fieldsCf)
        dataCF.(fieldsCf{ll})=dataCF.(fieldsCf{ll})(:,indTs);
    end

    %% Plot

    close all

    disp('Plotting widths ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);
    momentsSpecBasic.asl=HCRrange2asl(momentsSpecBasic.range,momentsSpecBasic.elevation,momentsSpecBasic.altitude);

    plotWidths(momentsTime,momentsSpecBasic,momentsSpecNoNoise,momentsSpecSmooth,momentsSpecCorrected,dataCF,figdir,project,showPlot);
    
end
