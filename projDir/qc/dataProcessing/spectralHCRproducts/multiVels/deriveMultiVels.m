% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.2';

% plotInds=0;
plotInds=(1:50:500);

outTime=0.1; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.1; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirTS(1:end-11),'figsMultiVel/cases/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/dualParticles_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

% For FIR filter
Fnorm=10/(500);           % Normalized frequency
firFilt=designfilt("lowpassfir",FilterOrder=70,CutoffFrequency=Fnorm);
filtShift=mean(grpdelay(firFilt));

for aa=1:length(caseStart)
    tic

    plotTimeAll=[];
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Moments
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.DBZ=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_CORR=[];
    dataCF.VEL_MASKED=[];
    % dataCF.eastward_velocity=[];
    % dataCF.northward_velocity=[];
    % dataCF.vertical_velocity=[];
    % dataCF.azimuth=[];
        
    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);

    % Find velocity correction for vel_raw to vel

    % xCorr=sind(dataCF.azimuth).*cosd(dataCF.elevation).*dataCF.eastward_velocity;
    % yCorr=cosd(dataCF.azimuth).*cosd(dataCF.elevation).*dataCF.northward_velocity;
    % zCorr=sind(dataCF.elevation).*dataCF.vertical_velocity;

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

        majorVelOne=single(nan(size(data.range,1),beamNum,1));
        minorVelOne=single(nan(size(data.range,1),beamNum,1));
        lowShoulderVelOne=single(nan(size(data.range,1),beamNum,1));
        highShoulderVelOne=single(nan(size(data.range,1),beamNum,1));

        peakVelsOne=single(nan(size(data.range,1),beamNum,3));
        peakPowsOne=single(nan(size(data.range,1),beamNum,3));

        % majorDbzOne=single(nan(size(data.range,1),beamNum,1));
        % minorDbzOne=single(nan(size(data.range,1),beamNum,1));
        % lowShoulderPowOne=single(nan(size(data.range,1),beamNum,1));
        % highShoulderPowOne=single(nan(size(data.range,1),beamNum,1));

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

            %% Spectral moments
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

            % This step removes the noise and also de-aliases
            [powerRMnoiseDBsmooth,specVelAdj,peakVels,peakPows]=noisePeaks(specPowerDB.V,momentsTimeOne.velRawDeAliased(:,ii),dataThis,firFilt,filtShift,plotTime);
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

            % Remove aircraft motion
            specVelAdj=specVelAdj+velBiasCorrection(:,cfInd);
            specVelRMnoise=specVelRMnoise+velBiasCorrection(:,cfInd);
            peakVels=peakVels+velBiasCorrection(:,cfInd);

            % Velocities at peaks            
            peakVelsOne(:,ii,:)=peakVels;
            peakPowsOne(:,ii,:)=peakPows;
           
            %% Find shoulders
            lowShoulderVelOne(:,ii)=specVelRMnoise(:,1);
            

            flipSpec=fliplr(specVelRMnoise);
            [~,highInds]=max(~isnan(flipSpec),[],2);
            highIndsLin=sub2ind(size(specVelRMnoise),1:size(specVelRMnoise,1),highInds');
            highShoulderVelOne(:,ii)=flipSpec(highIndsLin);

            %% Multi DBZ
            % % DBM
            % powerLinVmajor=10.^(majorDbzOne(:,ii,:)./10);
            % powerLinVminor=10.^(minorDbzOne(:,ii,:)./10);
            % 
            % % SNR
            % noiseLinV=10.^(dataThis.noise_v./10);
            % snrLinLowMajor=(powerLinVmajor-noiseLinV)./noiseLinV;
            % snrLinLowMajor(snrLinLowMajor<0)=nan;
            % snrDBmajor=10*log10(snrLinLowMajor);
            % 
            % snrLinLowMinor=(powerLinVminor-noiseLinV)./noiseLinV;
            % snrLinLowMinor(snrLinLowMinor<0)=nan;
            % snrDBminor=10*log10(snrLinLowMinor);
            % 
            % % DBZ
            % dataThis.range(dataThis.range<0)=nan;
            % majorDbzOne(:,ii,:)=snrDBmajor+20*log10(dataThis.range./1000)+data.dbz1km_v;
            % minorDbzOne(:,ii,:)=snrDBminor+20*log10(dataThis.range./1000)+data.dbz1km_v;
            % 
            % % Shoulders
            % % Linear
            % lowShoulderPowLin=10.^(lowShoulderPowOne./10);
            % highShoulderPowLin=10.^(highShoulderPowOne./10);
            % 
            % % SNR
            % snrLinLow=(lowShoulderPowLin-noiseLinV)./noiseLinV;
            % snrLinLow(snrLinLow<0)=nan;
            % snrLow=10*log10(snrLinLow);
            % snrLinHigh=(highShoulderPowLin-noiseLinV)./noiseLinV;
            % snrLinHigh(snrLinHigh<0)=nan;
            % snrHigh=10*log10(snrLinHigh);
            % 
            % % DBZ
            % data.range(data.range<0)=nan;
            % lowShoulderDbzOne=snrLow+20*log10(data.range./1000)+data.dbz1km_v;
            % highShoulderDbzOne=snrHigh+20*log10(data.range./1000)+data.dbz1km_v;
        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            majorVel=majorVelOne;
            minorVel=minorVelOne;
            lowShoulderVel=lowShoulderVelOne;
            highShoulderVel=highShoulderVelOne;
            peakVelsAll=peakVelsOne;
            peakPowsAll=peakPowsOne;

            % majorDbz=majorDbzOne;
            % minorDbz=minorDbzOne;
            % lowShoulderDbz=lowShoulderDbzOne;
            % highShoulderDbz=highShoulderDbzOne;
        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            lowShoulderVel=cat(2,lowShoulderVel,lowShoulderVelOne);
            highShoulderVel=cat(2,highShoulderVel,highShoulderVelOne);
            peakVelsAll=cat(2,peakVelsAll,peakVelsOne);
            peakPowsAll=cat(2,peakPowsAll,peakPowsOne);
            % lowShoulderDbz=cat(2,lowShoulderDbz,lowShoulderDbzOne);
            % highShoulderDbz=cat(2,highShoulderDbz,highShoulderDbzOne);

            % Major
            checkDims=size(majorVel,3)-size(majorVelOne,3);
            if checkDims<0
                majorVel=padarray(majorVel,[0,0,abs(checkDims)],nan,'post');
                % majorDbz=padarray(majorDbz,[0,0,abs(checkDims)],nan,'post');
            end

            momVelTest=nan(size(majorVelOne,1),size(majorVelOne,2),max(size(majorVel,3)));
            dim3=size(majorVelOne,3);
            momVelTest(:,:,1:dim3)=majorVelOne;     
            majorVel=cat(2,majorVel,momVelTest);

            % momPowTest=nan(size(majorVelOne,1),size(majorVelOne,2),max(size(majorVel,3)));
            % momPowTest(:,:,1:dim3)=majorDbzOne;     
            % majorDbz=cat(2,majorDbz,momPowTest);

            % Minor
            checkDims=size(minorVel,3)-size(minorVelOne,3);
            if checkDims<0
                minorVel=padarray(minorVel,[0,0,abs(checkDims)],nan,'post');
                % minorDbz=padarray(minorDbz,[0,0,abs(checkDims)],nan,'post');
            end

            momVelTest=nan(size(minorVelOne,1),size(minorVelOne,2),max(size(minorVel,3)));            
            dim3=size(minorVelOne,3);
            momVelTest(:,:,1:dim3)=minorVelOne;     
            minorVel=cat(2,minorVel,momVelTest);

            % momPowTest=nan(size(minorVelOne,1),size(minorVelOne,2),max(size(minorVel,3)));
            % momPowTest(:,:,1:dim3)=minorDbzOne;     
            % minorDbz=cat(2,minorDbz,momPowTest);
        end
       
    end
    
    %% Reverse up pointing direction
    shoulderLowVel=nan(size(lowShoulderVel));
    shoulderLowVel(:,momentsTime.elevation<=0)=lowShoulderVel(:,momentsTime.elevation<=0);
    shoulderLowVel(:,momentsTime.elevation>0)=-highShoulderVel(:,momentsTime.elevation>0);
    shoulderHighVel=nan(size(highShoulderVel));
    shoulderHighVel(:,momentsTime.elevation<=0)=highShoulderVel(:,momentsTime.elevation<=0);
    shoulderHighVel(:,momentsTime.elevation>0)=-lowShoulderVel(:,momentsTime.elevation>0);

    % shoulderLowDbz=nan(size(lowShoulderDbz));
    % shoulderLowDbz(:,momentsTime.elevation<=0)=lowShoulderDbz(:,momentsTime.elevation<=0);
    % shoulderLowDbz(:,momentsTime.elevation>0)=highShoulderDbz(:,momentsTime.elevation>0);
    % shoulderHighDbz=nan(size(highShoulderDbz));
    % shoulderHighDbz(:,momentsTime.elevation<=0)=highShoulderDbz(:,momentsTime.elevation<=0);
    % shoulderHighDbz(:,momentsTime.elevation>0)=lowShoulderDbz(:,momentsTime.elevation>0);

    % majorVel(:,momentsTime.elevation>0,:)=-majorVel(:,momentsTime.elevation>0,:);
    % minorVel(:,momentsTime.elevation>0,:)=-minorVel(:,momentsTime.elevation>0,:);

    peakVelsAll(:,momentsTime.elevation>0,:)=-peakVelsAll(:,momentsTime.elevation>0,:);

    % peakLowVel=nan(size(lowShoulderVel));
    % peakLowVel(:,momentsTime.elevation<=0)=peakVelsAll(:,momentsTime.elevation<=0,1);
    % peakLowVel(:,momentsTime.elevation>0)=-peakVelsAll(:,momentsTime.elevation>0,2);
    % peakHighVel=nan(size(highShoulderVel));
    % peakHighVel(:,momentsTime.elevation<=0)=peakVelsAll(:,momentsTime.elevation<=0,2);
    % peakHighVel(:,momentsTime.elevation>0)=-peakVelsAll(:,momentsTime.elevation>0,1);

    dataCF.VEL_MASKED(:,dataCF.elevation>0)=-dataCF.VEL_MASKED(:,dataCF.elevation>0);

    %% Sort vels into layers

    disp('Sorting layers ...')
    [velLayers,dbzLayers]=sortMultiVelLayers(peakVelsAll,peakPowsAll);
    
    %% Take time

    eSecs=toc;

    eData=momentsTime.time(end)-momentsTime.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Plot

    close all

    disp('Plotting velocities ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);

    plotMultiVels(momentsTime,dataCF,shoulderLowVel,shoulderHighVel,velLayers(:,:,1),velLayers(:,:,2),figdir,project,showPlot,plotTimeAll);
    %plotMultiRefs(momentsTime,shoulderLowDbz,shoulderHighDbz,dbzLayers,figdir,project,showPlot);

end
