% Analyze HCR time series

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.1';

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirTS,'figsAirVel/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airVel_',project,'.txt'];

freqStr=strfind(freqData,'hz');
outFreq=str2num(freqData(1:freqStr-1)); % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=7.5;
ylimLower=-0.1;

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

    dataTS=[];
    dataTS.IVc=[];
    dataTS.QVc=[];

    dataTS=readHCRts(fileListTS,dataTS,startTime,endTime);

    %% Load moments data
    disp("Getting moments data ...");

    fileListMoments=makeFileList(dataDirCF,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    dataCF=[];
    dataCF.DBZ=[];
    dataCF.VEL_RAW=[];
    dataCF.VEL_CORR=[];
    dataCF.VEL_MASKED=[];

    dataCF=read_HCR(fileListMoments,dataCF,startTime,endTime);

    dataCF.VEL_CORR(isnan(dataCF.VEL_MASKED))=nan;

    % Find velocity correction for vel_raw to vel

    velCorrection=dataCF.VEL_CORR-dataCF.VEL_RAW;
    velCorrection(abs(velCorrection)>5)=nan;
    velCorrection=median(velCorrection,1,'omitnan');

    %% Get correct samples

    beamInds=[];

    for ii=1:length(dataCF.time)
        [thisMin,thisInd]=min(abs(dataTS.time-dataCF.time(ii)));
        beamInds=[beamInds,thisInd];
    end

    beamInds=unique(beamInds,'stable');

    if length(beamInds)~=length(dataCF.time)
        error('Time series data missing')
    end

    sampleNum=floor((beamInds(3)-beamInds(2))./2)*2;

    startInds=beamInds-floor(sampleNum/2)+1;
    endInds=beamInds+floor(sampleNum/2);

    if startInds(1)<1
        startInds(1)=[];
        endInds(1)=[];
        beamInds(1)=[];
    end
    if endInds(end)>size(dataTS.IVc,2)
        startInds(end)=[];
        endInds(end)=[];
        beamInds(end)=[];
    end

    beamNum=length(startInds);

    %% Prepare processing
    % Nyquist velocity
    nyq=dataTS.lambda./(4*dataTS.prt(1));
    
    % Window
    win=window(@hamming,sampleNum);  % Default window is Hamming
    winWeight=sampleNum/sum(win);
    winNorm=win*winWeight;

    % Initialize output
    velDeAliasedSDall=nan(size(dataTS.range,1),beamNum);
    traceReflAll=nan(size(dataTS.range,1),beamNum);
    velAirAll=nan(size(dataTS.range,1),beamNum);
    velSmallerAll=nan(size(dataTS.range,1),beamNum);
    velLargerAll=nan(size(dataTS.range,1),beamNum);
    velMaxAll=nan(size(dataTS.range,1),beamNum);

    %% Loop through beams
    for ii=1:beamNum

        % Check if all nan       
        velDeAliased=dataCF.VEL_MASKED(:,ii);
        
        if sum(~isnan(velDeAliased))==0
            continue
        end

        %% Process spectra

        disp(['Beam ',num2str(ii),' of ',num2str(beamNum)]);

        % Start and end ind
        startInd=startInds(ii);
        endInd=endInds(ii);

        cIQv=winNorm'.*(dataTS.IVc(:,startInd:endInd)+i*dataTS.QVc(:,startInd:endInd))./sqrt(sampleNum);

        %% FFT and spectra

        fftIQ=fft(cIQv,[],2);

        powerRealIn=real(fftIQ).^2;
        powerImagIn=imag(fftIQ).^2;
        powerSignal=powerRealIn+powerImagIn;

        powerShifted=fftshift(powerSignal,2);

        % Reverse to get pointing direction consistent
        powerShifted=fliplr(powerShifted);

        powerSpec=10*log10(powerShifted);

        %% Adjust spectra so max is in the middle

        prtThis=mode(dataTS.prt(startInd:endInd));

        [powerAdj,phaseAdj]=adjSpecBounds(powerSpec,velDeAliased,sampleNum);

        %% Air velocity

        [velAir,velSmaller,velLarger,velMax,traceRefl]=getAirVel(powerAdj,phaseAdj,dataCF.elevation(ii),sampleNum,dataTS.lambda,prtThis,dataTS.range,dataTS.dbz1km_v,dataTS.noise_v);

        %% Velocity in spectral domain

        specLin=10.^(powerAdj./10);

        sumSpecLin=sum(specLin,2,'omitnan');
        sumSpecPhase=sum(specLin.*phaseAdj,2,'omitnan');

        meanK=sumSpecPhase./sumSpecLin;
        velSpec=dataTS.lambda/(4*pi.*prtThis).*meanK;

        %% Correct velocity for aircraft motion and bias

        velSpec=velSpec+velCorrection(ii);
        velAir=velAir+velCorrection(ii);
        velMax=velMax+velCorrection(ii);
        velLarger=velLarger+velCorrection(ii);
        velSmaller=velSmaller+velCorrection(ii);

        %% Find de-alias mask for spectral data

        deAliasDiff=velSpec-velDeAliased;

        deAliasMask=zeros(size(deAliasDiff));
        checkFold=[2,4,6];

        for jj=1:3
            deAliasMask(deAliasDiff>checkFold(jj)*nyq-3)=checkFold(jj)*nyq;
            deAliasMask(deAliasDiff<-(checkFold(jj)*nyq-3))=-checkFold(jj)*nyq;
        end

        %% Add de-aliasing and add to matrix

        if dataCF.elevation(ii)>0
            velDeAliasedSDall(:,ii)=-(velSpec-deAliasMask);
            velAirAll(:,ii)=-(velAir-deAliasMask);
            velSmallerAll(:,ii)=-(velSmaller-deAliasMask);
            velLargerAll(:,ii)=-(velLarger-deAliasMask);
            velMaxAll(:,ii)=-(velMax-deAliasMask);
        else
            velDeAliasedSDall(:,ii)=velSpec-deAliasMask;
            velAirAll(:,ii)=velAir-deAliasMask;
            velSmallerAll(:,ii)=velSmaller-deAliasMask;
            velLargerAll(:,ii)=velLarger-deAliasMask;
            velMaxAll(:,ii)=velMax-deAliasMask;
        end

        traceReflAll(:,ii)=traceRefl;

    end

    %% Plot vel
    close all

    dataCF.VEL_MASKED(:,dataCF.elevation>0)=-dataCF.VEL_MASKED(:,dataCF.elevation>0);

    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colM=colormap(velCols);
    colormap(colM);

    s1=subplot(3,1,1);
    surf(dataCF.time,dataCF.asl./1000,dataCF.VEL_MASKED,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity time domain de-aliased (m s^{-1})')

    s2=subplot(3,1,2);
    surf(dataCF.time,dataCF.asl./1000,velDeAliasedSDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity spectral domain de-aliased (m s^{-1})')

    s3=subplot(3,1,3);
    surf(dataCF.time,dataCF.asl./1000,velAirAll,'edgecolor','none');
    s3.Colormap=colM;
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Air velocity de-aliased (m s^{-1})')

    linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_vel_',datestr(dataCF.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(dataCF.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

    %% Plot air vel and reflectivities

    f1 = figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');

    s1=subplot(2,1,1);
    surf(dataCF.time,dataCF.asl./1000,dataCF.DBZ,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-40 20]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Reflectivity time domain (dBZ)');

    s2=subplot(2,1,2);
    surf(dataCF.time,dataCF.asl./1000,traceReflAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-40 20]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Tracer reflectivity spectral domain (dBZ)');

    linkaxes([s1 s2],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_dbz_',datestr(dataCF.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(dataCF.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

    %% Plot max, smaller, and larger vel

    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap(colM);

    s1=subplot(3,1,1);
    surf(dataCF.time,dataCF.asl./1000,velMaxAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity max peak (m s^{-1})')

    s2=subplot(3,1,2);
    surf(dataCF.time,dataCF.asl./1000,velLargerAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity larger vel peak (m s^{-1})')

    s3=subplot(3,1,3);
    surf(dataCF.time,dataCF.asl./1000,velSmallerAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([dataCF.time(1),dataCF.time(end)]);
    colorbar
    grid on
    box on
    title('Velocity smaller vel peak (m s^{-1})')

    linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velMaxLargerSmaller_',datestr(dataCF.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(dataCF.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end