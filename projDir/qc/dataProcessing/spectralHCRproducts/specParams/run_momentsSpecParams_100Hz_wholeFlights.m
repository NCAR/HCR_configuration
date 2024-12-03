% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc2';
freqData='10hz';
qcVersion='v2.0';
whichModel='era5';

saveData=1;

plotYes=1;
plotInds=0;
%plotInds=(1:50:500);

outTime=0.01; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.01; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'specMomentsParams/wholeFlights/'];

showPlot='off';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList = table2array(readtable(infile));

for aa=9:size(caseList,1)
    tic
    
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    plotTimeAll=[];

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
    fileListTS=makeFileList(dataDirTS,startTime,endTime,'20YYMMDDxhhmmss',1);

    % Check if first file is good
    firstGood=0;
    while firstGood==0
        data=[];
        data.IVc=[];
        data.QVc=[];
        try
            dataTest=read_TsArchive_iwrf_bulk(fileListTS{1},data);
            firstGood=1;
        catch
            fileListTS=fileListTS(2:end);
            warning('Removing first file.')
        end
    end

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

            try
                data=read_TsArchive_iwrf_bulk(fileListTS{bb},data);
            catch
                warning('Could not read file. Moving on.')
                continue
            end

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
            
            if minTimeDiff>0.1 | all(isnan(dataCF.VEL_MASKED(:,cfInd)))
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

    %% Combine to 10 hz

    disp('Averaging to 10 Hz ...');

    momentsSpecParams=[];
    momentsSpecParams.velRaw=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.width=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.skew=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.kurt=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.lrwidth=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.lslope=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.rslope=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.level=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.revel=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.lpvel=nan(size(dataCF.VEL_MASKED));
    momentsSpecParams.rpvel=nan(size(dataCF.VEL_MASKED));

    dataMask=double(~isnan(momentsSpec.skew));
    numNan=3; % Remove averages with fewer pixels

    fieldsTs=fieldnames(momentsSpecParams);
    for kk=1:length(dataCF.time)
        [minTD,matchTime]=min(abs(etime(datevec(momentsSpec.time),datevec(dataCF.time(kk)))));
        if minTD>0.05
            continue
        end
        tenInds=max([1,matchTime-4]):min([length(momentsSpec.time),matchTime+5]);
        dataMaskSum=sum(dataMask(:,tenInds),2);
        for ll=1:length(fieldsTs)
            momentsSpecParams.(fieldsTs{ll})(:,kk)=median(momentsSpec.(fieldsTs{ll})(:,tenInds),2,'omitnan');
            momentsSpecParams.(fieldsTs{ll})(dataMaskSum<numNan,kk)=nan;
        end
    end

    momentsSpecParams.time=dataCF.time;
    momentsSpecParams.asl=dataCF.asl;
    momentsSpecParams.range=dataCF.range;
    momentsSpecParams.elevation=dataCF.elevation;
    momentsSpecParams.longitude=dataCF.longitude;
    momentsSpecParams.latitude=dataCF.latitude;
    momentsSpecParams.altitude=dataCF.altitude;

    % Remove peak speckels
    momentsSpecParams.lpvel(isnan(momentsSpecParams.rpvel))=nan;
    momentsSpecParams.rpvel(isnan(momentsSpecParams.lpvel))=nan;

    rmRegs=50; % Minimum number of contiguous pixels
    pvelMask=~isnan(momentsSpecParams.lpvel);
    pvelMask=bwareaopen(pvelMask,rmRegs);

    momentsSpecParams.lpvel(pvelMask==0)=nan;
    momentsSpecParams.rpvel(pvelMask==0)=nan;

     %% Take time

    eSecs=toc;

    eData=momentsSpec.time(end)-momentsSpec.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);


    %% Plot in 20 minute increments

    if plotYes
        disp('Plotting ...');
        startPlot=startTime;

        fieldsPlot=fieldnames(momentsSpecParams);

        while startPlot<endTime

            endPlot=startPlot+minutes(20);
            indsTest=find(momentsSpecParams.time>=startPlot & momentsSpecParams.time<=endPlot);
            if length(indsTest)==0
                startPlot=endPlot;
                continue
            end
            
            newInds=indsTest(1):round(length(indsTest)/2000):indsTest(end);

            % Resample for plotting
            for ll=1:length(fieldsPlot)
                momentsSpecParamsPlot.(fieldsPlot{ll})=momentsSpecParams.(fieldsPlot{ll})(:,newInds);
            end

            if sum(sum(~isnan(momentsSpecParamsPlot.velRaw)))>300

                close all
                plotMomentsSpecParams(momentsSpecParamsPlot,figdir,project,showPlot);
            end
            startPlot=endPlot;
        end
    end
    
    %% Save
    if saveData
        disp('Saving moments and spectral parameters ...')
        momentsSpecParams.vel=momentsSpecParams.velRaw;
        momentsSpecParams=rmfield(momentsSpecParams,'velRaw');

        save([outdir,whichModel,'.momentsSpecParams.',datestr(momentsSpecParams.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(momentsSpecParams.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'momentsSpecParams','-v7.3');
    end
end
