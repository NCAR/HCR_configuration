% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF=[];
freqData=[];
qcVersion=[];

outTime=0.1; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.02; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDirTS(1:end-6),'figsMultiVel/cases/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/dualParticles_',project,'.txt'];

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

    if isempty(fileListTS)
        warning('No data found.');
        continue
    end

    for bb=1:length(fileListTS)

        if bb==1
            data=[];
            data.IVc=[];
            data.QVc=[];
            data.eastward_velocity=[];
            data.northward_velocity=[];
            data.vertical_velocity=[];
            data.azimuth_vc=[];

            vars=fieldnames(data);

            %% Load data TS
            disp(['Loading time series file ',num2str(bb),' of ' num2str(length(fileListTS)),' ...']);
            data=readHCRts(fileListTS(1),data,startTime-seconds(outTime),endTime+seconds(outTime),0);

            % Find available times
            timeTest=data.time';
            timeTest(data.time<startTime | data.time>endTime)=[];
            TTdata=timetable(timeTest,ones(size(timeTest)));
            synchTT=retime(TTdata,'regular','sum','TimeStep',seconds(outTime));
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
        momentsTimeOne.eastward_velocity=nan(1,beamNum);
        momentsTimeOne.northward_velocity=nan(1,beamNum);
        momentsTimeOne.vertical_velocity=nan(1,beamNum);
        momentsTimeOne.azimuth_vc=nan(1,beamNum);
        momentsTimeOne.time=goodTimes(1:beamNum)';

        majorVelOne=single(nan(size(data.range,1),beamNum,1));
        majorPowOne=single(nan(size(data.range,1),beamNum,1));
        minorVelOne=single(nan(size(data.range,1),beamNum,1));
        minorPowOne=single(nan(size(data.range,1),beamNum,1));

        lowShoulderVelOne=single(nan(size(data.range,1),beamNum,1));
        highShoulderVelOne=single(nan(size(data.range,1),beamNum,1));
        lowShoulderPowOne=single(nan(size(data.range,1),beamNum,1));
        highShoulderPowOne=single(nan(size(data.range,1),beamNum,1));

        % Loop through beams
        for ii=1:beamNum

            %disp(datestr(goodTimes(ii),'yyyymmdd_HHMMSS.FFF'));
            
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
            momentsTimeOne.eastward_velocity(ii)=median(dataThis.eastward_velocity);
            momentsTimeOne.northward_velocity(ii)=median(dataThis.northward_velocity);
            momentsTimeOne.vertical_velocity(ii)=median(dataThis.vertical_velocity);
            momentsTimeOne.azimuth_vc(ii)=median(dataThis.azimuth_vc);
            momentsTimeOne.altitude(ii)=median(dataThis.altitude);

            %% Time moments
            momentsTimeOne=calcMomentsTime(cIQ,ii,momentsTimeOne,dataThis);

            % Censor on V (find missing and NScal)
            if max(momentsTimeOne.powerV(1:14,ii))<-80 | (momentsTimeOne.powerV(1,ii)>-96 & momentsTimeOne.powerV(1,ii)<-87)
                censorY=nan(size(momentsTimeOne.range,1),1);
            else
                % Censor on SNR and NCP
                censorY=momentsTimeOne.snr(:,ii)<0 & momentsTimeOne.ncp(:,ii)<0.2;
                censorY(isnan(momentsTimeOne.snr(:,ii)))=1;
                censorY=~censorY;
                censorY=double(censorY);
                censorY(censorY==0)=nan;
                censorY=movmedian(censorY,7,'includemissing');
                censorY=movmedian(censorY,7,'omitmissing');
            end
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
            [velPrev,prevCount,prevKeep,flipYes]=setUpPrev(finalRay,velPrev,prevCount,prevKeep,flipYes,momentsTimeOne.elevation(ii),1/outTime,defaultPrev);
           
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

            [majorVelOne,majorPowOne,minorVelOne,minorPowOne,lowShoulderVelOne,highShoulderVelOne,lowShoulderPowOne,highShoulderPowOne] ...
                =findMultiVels(powerRMnoiseDBsmooth,specVelRMnoise,specPowerDBadj,majorVelOne,majorPowOne,minorVelOne,minorPowOne,lowShoulderVelOne,highShoulderVelOne,lowShoulderPowOne,highShoulderPowOne,ii);

            %% Correct for aircraft motion
            lowShoulderVelOne(:,ii)=lowShoulderVelOne(:,ii)+single(xCorr+yCorr+zCorr);
            highShoulderVelOne(:,ii)=highShoulderVelOne(:,ii)+single(xCorr+yCorr+zCorr);

            %% Multi DBZ
            % DBM
            powerLinV=10.^(majorPowOne(:,ii,:)./10);
            
            % SNR
            noiseLinV=10.^(dataThis.noise_v./10);
            snrLinV=(powerLinV-noiseLinV)./noiseLinV;
            snrLinV(snrLinV<0)=nan;
            snrDB=10*log10(snrLinV);

            % DBZ
            dataThis.range(dataThis.range<0)=nan;
            majorPowOne(:,ii,:)=snrDB+20*log10(dataThis.range./1000)+data.dbz1km_v;
        end

        %% Add to output
        if bb==1
            momentsTime=momentsTimeOne;
            majorVel=majorVelOne;
            majorPow=majorPowOne;
            minorVel=minorVelOne;
            minorPow=minorPowOne;
            lowShoulderVel=lowShoulderVelOne;
            highShoulderVel=highShoulderVelOne;
            lowShoulderPow=lowShoulderPowOne;
            highShoulderPow=highShoulderPowOne;
        else
            dataFields=fields(momentsTime);

            for hh=1:length(dataFields)
                momentsTime.(dataFields{hh})=cat(2,momentsTime.(dataFields{hh}),momentsTimeOne.(dataFields{hh}));
            end

            lowShoulderVel=cat(2,lowShoulderVel,lowShoulderVelOne);
            highShoulderVel=cat(2,highShoulderVel,highShoulderVelOne);
            lowShoulderPow=cat(2,lowShoulderPow,lowShoulderPowOne);
            highShoulderPow=cat(2,highShoulderPowOne,highShoulderPowOne);

            % Major
            checkDims=size(majorVel,3)-size(majorVelOne,3);
            if checkDims<0
                majorVel=padarray(majorVel,[0,0,abs(checkDims)],nan,'post');
                majorPow=padarray(majorPow,[0,0,abs(checkDims)],nan,'post');
            end

            momVelTest=nan(size(majorVelOne,1),size(majorVelOne,2),max(size(majorVel,3)));
            momPowTest=nan(size(majorVelOne,1),size(majorVelOne,2),max(size(majorVel,3)));
            dim3=size(majorVelOne,3);
            momVelTest(:,:,1:dim3)=majorVelOne;     
            majorVel=cat(2,majorVel,momVelTest);
            momPowTest(:,:,1:dim3)=majorPowOne;     
            majorPow=cat(2,majorPow,momPowTest);

            % Minor
            checkDims=size(minorVel,3)-size(minorVelOne,3);
            if checkDims<0
                minorVel=padarray(minorVel,[0,0,abs(checkDims)],nan,'post');
                minorPow=padarray(minorPow,[0,0,abs(checkDims)],nan,'post');
            end

            momVelTest=nan(size(minorVelOne,1),size(minorVelOne,2),max(size(minorVel,3)));
            momPowTest=nan(size(minorVelOne,1),size(minorVelOne,2),max(size(minorVel,3)));
            dim3=size(minorVelOne,3);
            momVelTest(:,:,1:dim3)=minorVelOne;     
            minorVel=cat(2,minorVel,momVelTest);
            momPowTest(:,:,1:dim3)=minorPowOne;     
            minorPow=cat(2,minorPow,momPowTest);
        end
    end
    
    %% Sort vels into layers
   [velLayers,powLayers]=sortVelLayers(majorVel,majorPow,minorVel,minorPow,momentsTime);

    %% Take time

    eSecs=toc;

    eData=momentsTime.time(end)-momentsTime.time(1);
    timePerMin=eSecs/60/minutes(eData);
    disp(['Total: ',num2str(eSecs/60),' minutes. Per data minute: ',num2str(timePerMin),' minutes.']);

    %% Plot

    close all

    disp('Plotting velocities ...');

    momentsTime.asl=HCRrange2asl(momentsTime.range,momentsTime.elevation,momentsTime.altitude);

    plotMultiVels(momentsTime,lowShoulderVel,highShoulderVel,figdir,project,showPlot);
    %plotReflectivities(momentsDbzDual,momentsTime,figdir,project,showPlot);

end
