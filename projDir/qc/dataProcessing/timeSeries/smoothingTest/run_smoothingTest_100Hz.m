% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc3';
freqData='10hz';
qcVersion='v3.2';

%plotInds=0;
plotInds=(1:50:500);

outTime=0.01; % Desired output time resolution in seconds. Must be less than or equal to one second.
sampleTime=0.01; % Length of sample in seconds.

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'spectralMoments/smoothingTest/100hz/'];

showPlot='on';

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/spectralSmoothing_',project,'.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,caseList.Var6);
caseEnd=datetime(caseList.Var7,caseList.Var8,caseList.Var9, ...
    caseList.Var10,caseList.Var11,caseList.Var12);

errAll=[];
residAll=[];
velAircAll=[];

errCases=cell(length(caseStart),1);
velAircCases=cell(length(caseStart),1);

for aa=1:length(caseStart)
    tic

    err=[];
    resid=[];
    velAirc=[];

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
            widthSquares=momentsTimeOne.width(:,ii).^2-widthCorrDelta(:,cfInd).^2;
            widthSquares(widthSquares<0.1)=0.01;
            momentsTimeOne.widthCorr(:,ii)=sqrt(widthSquares);

            %% Spectra
            [specPowerLin,specPowerDB]=getSpectra(cIQ);
          
             % Censor
            specPowerLin.V(isnan(momentsTimeOne.vel(:,ii)),:)=nan;
            specPowerDB.V(isnan(momentsTimeOne.vel(:,ii)),:)=nan;

            % % Censor on SNR
            % snrThresh=20;
            % specPowerDB.V(dataCF.SNR(:,cfInd)<snrThresh,:)=nan;

            %% Remove noise, find edge points and peaks
            if ismember(ii,plotInds)
                plotTime=momentsTimeOne.time(ii);
                plotTimeAll=cat(1,plotTimeAll,plotTime);
            else
                plotTime=[];
            end

            % This step removes the noise, de-aliases, (and corrects for
            % spectral broadening)
            [err,resid,velAirc]=noisePeaks_smoothingTest(specPowerDB.V, ...
                momentsTimeOne.velRawDeAliased(:,ii),dataThis,widthCorrDelta(:,cfInd),velTestWind(:,cfInd),err,resid,velAirc,figdir,plotTime);

            %velAirc=cat(2,velAirc,repmat(velHeadwind(:,cfInd),1,size(err,2)-size(velAirc,2)));
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

    errMean=mean(err,2);
    errStd=std(err,0,2);

    residStd=std(resid);

    [minErr,minInd]=min(errMean);
    [minStd,minStdInd]=min(errStd);

    numZero=2:250;

    bestZero=numZero(minInd);

    [minErrAll,minIndAll]=min(err,[],1);
    bestZeroAll=numZero(minIndAll);

    H=histcounts(bestZeroAll,numZero-0.5);

    [maxH,maxIndH]=max(H);
    
    %% Plot
    close all
    
    f1 = figure('Position',[200 500 700 850],'DefaultAxesFontSize',12,'renderer','painters');

    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
    s1=nexttile(1);

    hold on
    l1=plot(numZero,errMean,'-b','LineWidth',2);
    l2=plot(numZero,errMean+errStd,'-k','LineWidth',0.7);
    plot(numZero,errMean-errStd,'-k','LineWidth',0.7);
    ylim([4.5,7]);

    scatter(bestZero,minErr,60,'filled','red')
    scatter(numZero(maxIndH),errMean(maxIndH),60,'filled','green')

    text(50,6.8,['Number of spectra: ',num2str(size(err,2)/2),' x2'],'fontsize',12)
    text(50,6.6,['Error without smoothing: ',num2str(errMean(end),3)],'fontsize',12)
    text(50,6.4,['Error at peak of ',num2str(numZero(maxIndH)),' non-zeros: ',num2str(errMean(maxIndH))],'fontsize',12)
    text(50,6.2,['Minimum error at ',num2str(bestZero),' non-zeros: ',num2str(minErr)],'fontsize',12)
    text(50,6.0,['Minimum standard deviation at ',num2str(numZero(minStdInd)),' non-zeros: ',num2str(minStd)],'fontsize',12)

    text(50,4.9,['Standard deviation of noise: ',num2str(residStd)],'fontsize',12)
    text(50,4.7,['Mean aircraft speed: ',num2str(mean(velAirc)),' m s^{-1}'],'fontsize',12)

    ylabel('Root mean square error')

    yyaxis right
    l3=plot(numZero,errStd,'-c','LineWidth',2);
    ylim([0.2,0.45]);
    scatter(numZero(minStdInd),minStd,60,'filled','magenta')

    legend([l1,l2,l3],{'Mean','St. dev.','St. dev.'},'Location','southeast')

    grid on
    box on

    xlabel('Number of non-zeros')
    ylabel('Standard deviation')

    title('RMSE vs number of non-zeros')

    s2=nexttile(2);

    bar(numZero(1:end-1),H,1)

    xlabel('Number of non-zeros');
    title('Distribution of minima')

    grid on
    box on

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_smoothingAnalysis_everyOther_',datestr(momentsSpecBasic.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momentsSpecBasic.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

    %% Add to output
    errAll=cat(2,errAll,err);
    residAll=cat(2,residAll,resid);
    velAircAll=cat(2,velAircAll,velAirc);

    errCases{aa}=err;
    velAircCases{aa}=velAirc;

end

errMean=mean(errAll,2);
errStd=std(errAll,0,2);

residStd=std(residAll,0,2);

[minErr,minInd]=min(errMean);
[minStd,minStdInd]=min(errStd);

numZero=2:250;

bestZero=numZero(minInd);

[minErrAll,minIndAll]=min(errAll,[],1);
bestZeroAll=numZero(minIndAll);

H=histcounts(bestZeroAll,numZero-0.5);

[maxH,maxIndH]=max(H);

%% Plot
close all

f1 = figure('Position',[200 500 700 850],'DefaultAxesFontSize',12,'renderer','painters');

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

hold on
l1=plot(numZero,errMean,'-b','LineWidth',2);
l2=plot(numZero,errMean+errStd,'-k','LineWidth',0.7);
plot(numZero,errMean-errStd,'-k','LineWidth',0.7);
ylim([4.5,7]);

scatter(bestZero,minErr,60,'filled','red')
scatter(numZero(maxIndH),errMean(maxIndH),60,'filled','green')

text(50,6.8,['Number of spectra: ',num2str(size(errAll,2)/2),' x2'],'fontsize',12)
text(50,6.6,['Error without smoothing: ',num2str(errMean(end),3)],'fontsize',12)
text(50,6.4,['Error at peak of ',num2str(numZero(maxIndH)),' non-zeros: ',num2str(errMean(maxIndH))],'fontsize',12)
text(50,6.2,['Minimum error at ',num2str(bestZero),' non-zeros: ',num2str(minErr)],'fontsize',12)
text(50,6.0,['Minimum standard deviation at ',num2str(numZero(minStdInd)),' non-zeros: ',num2str(minStd)],'fontsize',12)

text(50,4.7,['Standard deviation of noise: ',num2str(residStd)],'fontsize',12)

ylabel('Root mean square error')

yyaxis right
l3=plot(numZero,errStd,'-c','LineWidth',2);
ylim([0.2,0.45]);
scatter(numZero(minStdInd),minStd,60,'filled','magenta')

legend([l1,l2,l3],{'Mean','St. dev.','St. dev.'},'Location','southeast')

grid on
box on

xlabel('Number of non-zeros')
ylabel('Standard deviation')

title('RMSE vs number of non-zeros')

s2=nexttile(2);

bar(numZero(1:end-1),H,1)

xlabel('Number of non-zeros');
title('Distribution of minima')

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_smoothingAnalysis_everyOther'],'-dpng','-r0');

%% Split up by aircraft speed

minSpeed=min(velAircAll);
maxSpeed=max(velAircAll);

startInt=floor(minSpeed/10)*10;
endInt=ceil(maxSpeed/10)*10;

edgesInt=startInt:10:endInt;
edgesHalf=edgesInt(1:end-1)+5;

peakNum=nan(length(edgesInt)-1,1);

for ii=1:length(edgesInt)-1

    speedInds=find(velAircAll>edgesInt(ii) & velAircAll<=edgesInt(ii+1));
    if isempty(speedInds)
        continue
    end
    errAllInt=errAll(:,speedInds);
    [minErrAllInt,minIndAllInt]=min(errAllInt,[],1);
    bestZeroAllInt=numZero(minIndAllInt);

    Hint=histcounts(bestZeroAllInt,numZero-0.5);

    [maxHint,maxIndHint]=max(Hint);

    peakNum(ii)=numZero(maxIndHint);
end

%% Scatter plot
close all

f1 = figure('Position',[200 500 700 500],'DefaultAxesFontSize',12,'renderer','painters');

scatter(edgesHalf,peakNum,'filled');
xlim([edgesInt(1),edgesInt(end)]);
ylim([min(peakNum)-1,max(peakNum)+1]);

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['Ideal smoothing vs aircraft speed: ',project]);

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_smoothingAnalysis_everyOther_aircraftSpeed'],'-dpng','-r0');

save([figdir,project,'_smoothingAnalysis_everyOther_aircraftSpeed.mat'],'edgesHalf','peakNum');

%% Plot cases

% Collect data for scatter plot
peakNumCases=nan(length(errCases),1);
velAircMeanC=nan(length(errCases),1);

for ii=1:length(errCases)
    errC=errCases{ii};
    [minErrC,minIndErrC]=min(errC,[],1);
    bestZeroErrC=numZero(minIndErrC);

    HC=histcounts(bestZeroErrC,numZero-0.5);

    [maxHC,maxIndHC]=max(HC);

    peakNumCases(ii)=numZero(maxIndHC);

    velAircC=velAircCases{ii};
    velAircMeanC(ii)=mean(velAircC);
end

f2 = figure('Position',[200 500 700 500],'DefaultAxesFontSize',12,'renderer','painters');

scatter(velAircMeanC,peakNumCases,'filled');
xlim([edgesInt(1),edgesInt(end)]);
ylim([min(peakNumCases)-1,max(peakNumCases)+1]);

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['Ideal smoothing vs aircraft speed for cases: ',project]);

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f2,[figdir,project,'_smoothingAnalysis_everyOther_aircraftSpeed_cases'],'-dpng','-r0');

save([figdir,project,'_smoothingAnalysis_everyOther_aircraftSpeed_cases.mat'],'velAircMeanC','peakNumCases');