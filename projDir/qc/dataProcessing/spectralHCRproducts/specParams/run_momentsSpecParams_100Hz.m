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
%plotInds=(1:500:5000);

outTime=0.01; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.01; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'specMomentsParams/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airMotion_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

for aa=11:length(caseStart)
    tic

    plotTimeAll=[];
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    %startTime=datetime(2021,5,29,19,11,33);
    endTime=caseEnd(aa);

    %% CfRadial Moments
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_MASKED=[];
    dataCF.eastward_velocity=[];
    dataCF.northward_velocity=[];
         
    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);
    dataCF.beamWidth=ncread(fileListMoments{1},'radar_beam_width_v');

    % Find width correction factor and filter value
    % Aircraft speed
    velTestWind=sqrt(dataCF.eastward_velocity.^2+dataCF.northward_velocity.^2);
    % Correction factor
    widthCorrDelta=abs(0.3.*velTestWind.*sin(deg2rad(dataCF.elevation)).*deg2rad(dataCF.beamWidth));
    widthCorrDelta=fillmissing(widthCorrDelta,'nearest',1);
    widthCorrDelta=repmat(widthCorrDelta,size(dataCF.range,1),1);
        
    % Filter value
    if sampleTime==0.1
        filterAt=round(0.00022396.*velTestWind.^2-0.10542.*velTestWind+18.132);
    elseif sampleTime==0.01
        filterAt=round(0.000126.*velTestWind.^2-0.0548.*velTestWind+10.7);
    else
        error('Sample time must be 0.1 or 0.01.')
    end
    
    filterAt=fillmissing(filterAt,'nearest');
    filterAt=round(movmean(filterAt,501));
    filterAt=modefilt(filterAt,[1,101]);

    filterAt=repmat(filterAt,size(dataCF.range,1),1);
    
    % Velocity bias term
    velBiasCorrection=dataCF.VEL_MASKED-dataCF.VEL_RAW;
        
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

            % Bias correction
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
        momentsSpecOne.velRaw=nan(size(data.range,1),beamNum);
        momentsSpecOne.width=nan(size(data.range,1),beamNum);
        momentsSpecOne.skew=nan(size(data.range,1),beamNum);
        momentsSpecOne.kurt=nan(size(data.range,1),beamNum);
        momentsSpecOne.lrwidth=nan(size(data.range,1),beamNum);
        momentsSpecOne.lslope=nan(size(data.range,1),beamNum);
        momentsSpecOne.rslope=nan(size(data.range,1),beamNum);
        momentsSpecOne.level=nan(size(data.range,1),beamNum);
        momentsSpecOne.revel=nan(size(data.range,1),beamNum);
        momentsSpecOne.lpvel=nan(size(data.range,1),beamNum);
        momentsSpecOne.rpvel=nan(size(data.range,1),beamNum);
        momentsSpecOne.time=goodTimes(1:beamNum)';

        noiseFloor=nan(size(data.range,1),beamNum);
        velForDeAliasOne=nan(size(data.range,1),beamNum);

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
                        
            %% Find time index in CF moments
            [minTimeDiff,cfInd]=min(abs(etime(datevec(dataCF.time),datevec(momentsSpecOne.time(ii)))));
            
            if minTimeDiff>0.1
                continue
            end

            %% Spectra
            [specPowerLin,specPowerDB]=getSpectra(cIQ);
          
             % Censor
            specPowerLin.V(isnan(dataCF.VEL_MASKED(:,cfInd)),:)=nan;
            specPowerDB.V(isnan(dataCF.VEL_MASKED(:,cfInd)),:)=nan;

            %% Remove noise, find edge points and peaks
            if ismember(ii,plotInds)
                plotTime=momentsSpecOne.time(ii);
                plotTimeAll=cat(1,plotTimeAll,plotTime);
            else
                
                plotTime=[];
            end

            % This step removes the noise, de-aliases, and corrects for
            % spectral broadening
            [powerOrig,powerOrigRMnoise,powerSmooth,powerSmoothCorr,specVelAdj,noiseFloor(:,ii),peaks1,peaks2]= ...
                noisePeaks_skewKurtSP(specPowerDB.V,dataThis,widthCorrDelta(:,cfInd),filterAt(:,cfInd), ...
                sampleTime,figdir,plotTime);
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerSmoothCorr))=nan;

            % Remove aircraft motion
            specVelRMnoise=specVelRMnoise+velBiasCorrection(:,cfInd);
            
            %% Spectral moments

            momentsSpecOne=calcMomentsSpec_higherMoments(powerSmoothCorr,specVelRMnoise,ii,momentsSpecOne);
            
             %% Spectral parameters

            momentsSpecOne=calcSpecParams(powerSmoothCorr,specVelRMnoise,peaks1,peaks2,noiseFloor(:,ii),ii,momentsSpecOne);
            
            %% Set up de-aliasing
            velForDeAliasOne(:,ii)=dataCF.VEL_MASKED(:,cfInd);
        end

        %% Add to output
        if bb==1
            momentsSpec=momentsSpecOne;
            velForDeAlias=velForDeAliasOne;
        else
            dataFields1=fields(momentsSpecOne);

            for hh=1:length(dataFields1)
                momentsSpec.(dataFields1{hh})=cat(2,momentsSpec.(dataFields1{hh}),momentsSpecOne.(dataFields1{hh}));
            end

            velForDeAlias=cat(2,velForDeAlias,velForDeAliasOne);
        end
       
    end

     %% De-alias

    velDeAliasCheck=velForDeAlias-momentsSpec.velRaw;

    deAliasMaskE=zeros(size(velDeAliasCheck));
    for jj=1:3
        deAliasMaskE(velDeAliasCheck>checkFold(jj)*toVel-5)=checkFold(jj)*toVel;
        deAliasMaskE(velDeAliasCheck<-(checkFold(jj)*toVel-5))=-checkFold(jj)*toVel;
    end

    momentsSpec.velRaw=momentsSpec.velRaw+deAliasMaskE;
    momentsSpec.level=momentsSpec.level+deAliasMaskE;
    momentsSpec.revel=momentsSpec.revel+deAliasMaskE;
    momentsSpec.lpvel=momentsSpec.lpvel+deAliasMaskE;
    momentsSpec.rpvel=momentsSpec.rpvel+deAliasMaskE;

    %% Take time

    eSecs=toc;

    eData=momentsSpec.time(end)-momentsSpec.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Combine to 10 hz

    momentsSC10=[];
    momentsSC10.velRaw=nan(size(dataCF.VEL_MASKED));
    momentsSC10.width=nan(size(dataCF.VEL_MASKED));
    momentsSC10.skew=nan(size(dataCF.VEL_MASKED));
    momentsSC10.kurt=nan(size(dataCF.VEL_MASKED));
    momentsSC10.lrwidth=nan(size(dataCF.VEL_MASKED));
    momentsSC10.lslope=nan(size(dataCF.VEL_MASKED));
    momentsSC10.rslope=nan(size(dataCF.VEL_MASKED));
    momentsSC10.level=nan(size(dataCF.VEL_MASKED));
    momentsSC10.revel=nan(size(dataCF.VEL_MASKED));
    momentsSC10.lpvel=nan(size(dataCF.VEL_MASKED));
    momentsSC10.rpvel=nan(size(dataCF.VEL_MASKED));

    dataMask=double(~isnan(momentsSpec.skew));
    numNan=3; % Remove averages with fewer pixels

    fieldsTs=fieldnames(momentsSC10);
    for kk=1:length(dataCF.time)
        [minTD,matchTime]=min(abs(etime(datevec(momentsSpec.time),datevec(dataCF.time(kk)))));
        if minTD>0.05
            continue
        end
        tenInds=max([1,matchTime-4]):min([length(momentsSpec.time),matchTime+5]);
        dataMaskSum=sum(dataMask(:,tenInds),2);
        for ll=1:length(fieldsTs)
            momentsSC10.(fieldsTs{ll})(:,kk)=median(momentsSpec.(fieldsTs{ll})(:,tenInds),2,'omitnan');
            momentsSC10.(fieldsTs{ll})(dataMaskSum<numNan,kk)=nan;
        end
    end

    momentsSC10.time=dataCF.time;
    momentsSC10.asl=dataCF.asl;
    momentsSC10.range=dataCF.range;
    momentsSC10.elevation=dataCF.elevation;
    momentsSC10.longitude=dataCF.longitude;
    momentsSC10.latitude=dataCF.latitude;
    momentsSC10.altitude=dataCF.altitude;

    % Remove peak speckels
    momentsSC10.lpvel(isnan(momentsSC10.rpvel))=nan;
    momentsSC10.rpvel(isnan(momentsSC10.lpvel))=nan;

    rmRegs=50; % Minimum number of contiguous pixels
    pvelMask=~isnan(momentsSC10.lpvel);
    pvelMask=bwareaopen(pvelMask,rmRegs);

    momentsSC10.lpvel(pvelMask==0)=nan;
    momentsSC10.rpvel(pvelMask==0)=nan;

    %% Plot

    close all

    disp('Plotting ...');

    plotMomentsSpecParams(momentsSC10,figdir,project,showPlot);
end
