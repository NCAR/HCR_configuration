% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc3';
freqData='10hz'; % !!!!!!!!! Must be equal or less than one second !!!!!!!!!!!!!
qcVersion='v3.2';

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
    tic

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    for bb=1:length(fileListTS)

        if bb==1
            data=[];
            data.IVc=[];
            data.QVc=[];
            data.IHx=[];
            data.QHx=[];
            data.eastward_velocity=[];
            data.northward_velocity=[];
            data.vertical_velocity=[];
            data.azimuth_vc=[];

            vars=fieldnames(data);

            %% Load data TS
            disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);
            data=readHCRts(fileListTS(1),data,startTime-seconds(timeSpan),endTime+seconds(timeSpan),0);

            % Find available times
            timeTest=data.time';
            timeTest(data.time<startTime | data.time>endTime)=[];
            TTdata=timetable(timeTest,ones(size(timeTest)));
            synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(timeSpan));
            goodTimes=synchTT.timeTest(synchTT.Var1>0);

            beamNum=length(goodTimes)-1;

            % Set up de-aliasing
            defaultPrev=data.lambda/(mode(data.prt)*4);

            velPrev=repmat(defaultPrev,size(data.range,1),1);
            prevCount=zeros(size(velPrev));
            prevKeep=nan(size(velPrev));
            flipYes=0;
        else
            %% Trimm first file
            data=trimFirstFile(data,endInd);

            goodTimes(1:end-1)=[];

            %% Load data TS
            disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);

            data=readHCRts_add(fileListTS{bb},data,vars);

            % Find available times
            timeTest=data.time';
            timeTest(data.time<goodTimes | data.time>endTime)=[];
        end

        TTdata=timetable(timeTest,ones(size(timeTest)));
        synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(timeSpan));
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
        momentsTimeOne.eastward_velocity=nan(1,beamNum);
        momentsTimeOne.northward_velocity=nan(1,beamNum);
        momentsTimeOne.vertical_velocity=nan(1,beamNum);
        momentsTimeOne.azimuth_vc=nan(1,beamNum);
        momentsTimeOne.time=goodTimes(1:beamNum)';

        momentsVelDualRawOne=single(nan(size(data.range,1),beamNum,1));

        % Loop through beams
        for ii=1:beamNum

            %disp(datestr(goodTimes(ii),'yyyymmdd_HHMMSS.FFF'));
            
            % Find start and end indices for beam
            [~,startInd]=min(abs(goodTimes(ii)-seconds(timeSpan/10)-data.time)); % Default: 20
            [~,endInd]=min(abs(goodTimes(ii)+seconds(timeSpan/10)-data.time));

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
            momentsTimeOne.range(:,ii)=dataThis.range;
            momentsTimeOne.elevation(ii)=median(dataThis.elevation);
            momentsTimeOne.eastward_velocity(ii)=median(dataThis.eastward_velocity);
            momentsTimeOne.northward_velocity(ii)=median(dataThis.northward_velocity);
            momentsTimeOne.vertical_velocity(ii)=median(dataThis.vertical_velocity);
            momentsTimeOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsTimeOne.altitude(ii)=median(dataThis.altitude);

            %% Time moments
            momentsTimeOne=calcMomentsTime(cIQ,ii,momentsTimeOne,dataThis);

            % Censor on SNR and NCP
            censorY=momentsTimeOne.snr(:,ii)<0 & momentsTimeOne.ncp(:,ii)<0.2;
            censorY(isnan(momentsTimeOne.snr(:,ii)))=1;
            censorY=~censorY;
            censorY=double(censorY);
            censorY(censorY==0)=nan;
            censorY=movmedian(censorY,7,'includemissing');
            censorY=movmedian(censorY,7,'omitmissing');
            momentsTimeOne.velRaw(isnan(censorY),ii)=nan;

            %% Correct time domain velocity folding

            if momentsTimeOne.elevation(ii)>0
                velRay=-momentsTimeOne.velRaw(:,ii);
            else
                velRay=momentsTimeOne.velRaw(:,ii);
            end
            velRay(1:12)=nan;

            finalRay=deAliasSingleRay(velRay,velPrev,defaultPrev,[],momentsTimeOne.time(ii));

            % Set up time consistency check

            [velPrev,prevCount,prevKeep,flipYes]=setUpPrev(finalRay,velPrev,prevCount,prevKeep,flipYes,momentsTimeOne.elevation(ii),outFreq,defaultPrev);

            
            %% Add to output
            if momentsTimeOne.elevation(ii)>0
                momentsTimeOne.velRawDeAliased(:,ii)=-finalRay;
            else
                momentsTimeOne.velRawDeAliased(:,ii)=finalRay;
            end

            %% Correct velocity for aircraft motion
            xCorr=sind(momentsTimeOne.azimuth_vc(ii)).*cosd(momentsTimeOne.elevation(ii)).*momentsTimeOne.eastward_velocity(ii);
            yCorr=cosd(momentsTimeOne.azimuth_vc(ii)).*cosd(momentsTimeOne.elevation(ii)).*momentsTimeOne.northward_velocity(ii);
            zCorr=sind(momentsTimeOne.elevation(ii)).*momentsTimeOne.vertical_velocity(ii);
            momentsTimeOne.vel(:,ii)=momentsTimeOne.velRawDeAliased(:,ii)+single(xCorr+yCorr+zCorr);

            %% Spectral moments
            [specPowerLin,specPowerDB]=getSpectra(cIQ);

            % Move peak of spectra to middle
            [specPowerDBadj,specVelAdj]=adjSpecBoundsV(specPowerDB.V,momentsTimeOne.velRawDeAliased(:,ii),dataThis);

            % De-alias
            % Mean vel
            noiseLinV=10.^(data.noise_v./10);
            x=specVelAdj;
            y=specPowerLin.V-noiseLinV;
            velMean=sum(y.*x,2,'omitnan')./sum(y,2,'omitnan');

            deAliasDiff=velMean-momentsTimeOne.velRawDeAliased(:,ii);

            deAliasMask=zeros(size(deAliasDiff));
            checkFold=[2,4,6];

            for jj=1:3
                deAliasMask(deAliasDiff>checkFold(jj)*defaultPrev-5)=checkFold(jj)*defaultPrev;
                deAliasMask(deAliasDiff<-(checkFold(jj)*defaultPrev-5))=-checkFold(jj)*defaultPrev;
            end

            specVelAdj=specVelAdj-deAliasMask;
          
            % Censor
            specPowerDBadj(isnan(momentsTimeOne.vel(:,ii)),:)=nan;
            specVelAdj(isnan(momentsTimeOne.vel(:,ii)),:)=nan;

            %% Remove noise
            [powerDBsmooth,powerRMnoiseDBsmooth]=rmNoiseSpec(specPowerDBadj);
            powerRMnoiseDB=specPowerDBadj;
            powerRMnoiseDB(isnan(powerRMnoiseDBsmooth))=nan;
            specVelRMnoise=specVelAdj;
            specVelRMnoise(isnan(powerRMnoiseDBsmooth))=nan;

            %% Find regions with dual particle species

            momentsVelDualRawOne=findDualParticles(powerRMnoiseDBsmooth,specVelRMnoise,specPowerDBadj,momentsVelDualRawOne,ii);

        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            momentsVelDualRaw=momentsVelDualRawOne;
        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            checkDims=size(momentsVelDualRaw,3)-size(momentsVelDualRawOne,3);
            if checkDims<0
                momentsVelDualRaw=padarray(momentsVelDualRaw,[0,0,abs(checkDims)],nan,'post');
            end

            momVelTest=nan(size(momentsVelDualRawOne,1),size(momentsVelDualRawOne,2),max(size(momentsVelDualRaw,3)));
            dim3=size(momentsVelDualRawOne,3);
            momVelTest(:,:,1:dim3)=momentsVelDualRawOne;     
            momentsVelDualRaw=cat(2,momentsVelDualRaw,momVelTest);
        end
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
