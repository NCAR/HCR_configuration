% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset, otrec
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

figdir=[dataDirCF(1:end-5),'specParams/'];

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
    %startTime=datetime(2021,5,29,19,11,11);
    endTime=caseEnd(aa);

    %% CfRadial Moments
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.DBZ=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_CORR=[];
    dataCF.VEL_MASKED=[];
    dataCF.LDR=[];
    dataCF.SNR=[];
    % dataCF.U=[];
    % dataCF.V=[];
    dataCF.eastward_velocity=[];
    dataCF.northward_velocity=[];
         
    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);
    dataCF.beamWidth=ncread(fileListMoments{1},'radar_beam_width_v');

    % Find width correction factor
    % Aircraft speed
    velTestWind=sqrt(dataCF.eastward_velocity.^2+dataCF.northward_velocity.^2);
    widthCorrDelta=abs(0.3.*velTestWind.*sin(deg2rad(dataCF.elevation)).*deg2rad(dataCF.beamWidth));
    widthCorrDelta=fillmissing(widthCorrDelta,'nearest',1);
    widthCorrDelta=repmat(widthCorrDelta,size(dataCF.range,1),1);  
    velTestWind=repmat(velTestWind,size(dataCF.range,1),1);

    % Or headwind
    % uHeadwind=dataCF.eastward_velocity-dataCF.U;
    % vHeadwind=dataCF.northward_velocity-dataCF.V;
    % velTestWind=sqrt(uHeadwind.^2+vHeadwind.^2);
    % widthCorrDelta=abs(0.3.*velTestWind.*sin(deg2rad(dataCF.elevation)).*deg2rad(dataCF.beamWidth));
    
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

        momentsSpecSmoothCorrOne=momentsTimeOne;
        momentsSpecSmoothCorrOne.lrwidth=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.lslope=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.rslope=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.level=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.revel=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.lpvel=nan(size(data.range,1),beamNum);
        momentsSpecSmoothCorrOne.rpvel=nan(size(data.range,1),beamNum);
        
        noiseFloor=nan(size(data.range,1),beamNum);

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
            %win=window(@tukeywin,sampleNum,0.5); % Tukey window
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

            momentsSpecSmoothCorrOne.range(:,ii)=dataThis.range;
            momentsSpecSmoothCorrOne.elevation(ii)=median(dataThis.elevation);
            momentsSpecSmoothCorrOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsSpecSmoothCorrOne.altitude(:,ii)=median(dataThis.altitude);

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
            widthSquares=momentsTimeOne.width(:,ii).^2-widthCorrDelta(:,cfInd).^2;
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

            % This step removes the noise, de-aliases, and corrects for
            % spectral broadening
            [powerOrig,powerOrigRMnoise,powerSmooth,powerSmoothCorr,specVelAdj,noiseFloor(:,ii),peaks1,peaks2]= ...
                noisePeaks_smoothCorr(specPowerDB.V,momentsTimeOne.velRawDeAliased(:,ii), ...
                dataThis,widthCorrDelta(:,cfInd),velTestWind(:,cfInd),sampleTime,figdir,plotTime);
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerSmoothCorr))=nan;

            % Remove aircraft motion
            specVelAdj=specVelAdj+velBiasCorrection(:,cfInd);
            specVelRMnoise=specVelRMnoise+velBiasCorrection(:,cfInd);
            
            %% Spectral moments

            momentsSpecSmoothCorrOne=calcMomentsSpec_higherMoments(powerSmoothCorr,specVelRMnoise,ii,momentsSpecSmoothCorrOne);

            %% Spectral parameters
            momentsSpecSmoothCorrOne=calcSpecParams(powerSmoothCorr,specVelRMnoise,peaks1,peaks2,noiseFloor(:,ii),ii,momentsSpecSmoothCorrOne);
            
        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            momentsSpecSmoothCorr=momentsSpecSmoothCorrOne;
        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            dataFields1=fields(momentsSpecSmoothCorrOne);

            for hh=1:length(dataFields1)
                momentsSpecSmoothCorr.(dataFields1{hh})=cat(2,momentsSpecSmoothCorr.(dataFields1{hh}),momentsSpecSmoothCorrOne.(dataFields1{hh}));
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
    end
    fieldsSpec=fieldnames(momentsSpecSmoothCorr);
    for ll=1:length(fieldsSpec)
        momentsSpecSmoothCorr.(fieldsSpec{ll})=momentsSpecSmoothCorr.(fieldsSpec{ll})(:,indTs);
    end

    dataCF=rmfield(dataCF,'beamWidth');

    fieldsCf=fieldnames(dataCF);
    for ll=1:length(fieldsCf)
        dataCF.(fieldsCf{ll})=dataCF.(fieldsCf{ll})(:,indTs);
    end

    %% Plot

    close all

    disp('Plotting ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);
    momentsSpecSmoothCorr.asl=HCRrange2asl(momentsSpecSmoothCorr.range,momentsSpecSmoothCorr.elevation,momentsSpecSmoothCorr.altitude);

    plotAllMoments(dataCF,momentsSpecSmoothCorr,figdir,project,showPlot);   
    plotSpecParams(momentsSpecSmoothCorr,figdir,project,showPlot);
    
end
