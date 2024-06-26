% Analyze HCR time series
clear all;
close all;

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.1';

dataDirTS=HCRdir(project,quality,qcVersion,freqData);
dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

figdir=[dataDirCF(1:end-5),'figsAirVel/spectral/'];

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

for aa=3:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Load data TS
    disp("Getting time series data ...");

    fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

    data=[];
    data.IVc=[];
    data.QVc=[];

    data=readHCRts(fileListTS,data,startTime,endTime);

    %% Calculate moments
    beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

    powerV=nan(size(data.range,1),beamNum);
    vel=nan(size(data.range,1),beamNum);
    width=nan(size(data.range,1),beamNum);
    dbz=nan(size(data.range,1),beamNum);
    snr=nan(size(data.range,1),beamNum);
    skew=nan(size(data.range,1),beamNum);
    kurt=nan(size(data.range,1),beamNum);

    timeBeams=[];

    startInd=1;
    endInd=1;
    ii=1;

    % Loop through beams
    while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

        % Find start and end indices for beam
        timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
        [minDiff,endInd]=min(abs(timeDiff+timeSpan));

        sampleNum=endInd-startInd+1;

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

        %prtThis=mean(prt(startInd:endInd));
        prt=mode(data.prt);

        %% FFT and spectra

        fftIQ=fft(cIQv,[],2);

        powerRealIn=real(fftIQ).^2;
        powerImagIn=imag(fftIQ).^2;
        powerSignal=powerRealIn+powerImagIn;

        powerShifted=fftshift(powerSignal,2);

        % Reverse to get pointing direction consistent
        specLinOrig=fliplr(powerShifted);

        % Find peak
        specLin=nan(size(specLinOrig));
        specVelVec=nan(size(specLinOrig));
        %specVelVecOrig=-pi:2*pi/(sampleNum-1):pi;
        specVelVecOrig=-3*pi:2*pi/(sampleNum):3*pi;
        specVelVecOrig=specVelVecOrig(1:end-1);

        for kk=1:size(specLinOrig,1)
            [~,maxInd]=max(specLinOrig(kk,:),[],'omitnan');
            maxInd=maxInd+sampleNum;
            sBsSpec=repmat(specLinOrig(kk,:),1,3);

            try
                specLin(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
                specVelVec(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
            catch
                specLin(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
                specVelVec(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
            end
        end

        %%%%%%%%%%%%%%%%%%

        % DBM
        powerLin=mean(specLin,2);
        powerV(:,ii)=10*log10(powerLin)-data.rx_gain_v;

        % SNR
        noiseLin=10.^(data.noise_v./10);
        snrLin=(powerLin-noiseLin)./noiseLin;
        snrLin(snrLin<0)=nan;
        snr(:,ii)=10*log10(snrLin);

        % DBZ
        data.range(data.range<0)=nan;
        dbz(:,ii)=snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

        % VEL
        vel(:,ii)=sum(specLin.*specVelVec,2,'omitnan')./sum(specLin,2,'omitnan');

        % WIDTH
        width(:,ii)=(sum(specLin.*(specVelVec-vel(:,ii)).^2,2,'omitnan')./sum(specLin,2,'omitnan')).^0.5;

        % SKEWNESS
        skew(:,ii)=sum(specLin.*(specVelVec-vel(:,ii)).^3,2,'omitnan')./(sum(specLin,2,'omitnan').*width(:,ii).^3);

        % KURTOSIS
        kurt(:,ii)=sum(specLin.*(specVelVec-vel(:,ii)).^4,2,'omitnan')./(sum(specLin,2,'omitnan').*width(:,ii).^4);

        timeBeams=[timeBeams;data.time(startInd)];

        if data.elevation(startInd)<0
            flipYes=1;
        else
            flipYes=0;
        end

        %% Plot
        if abs(etime(datevec(data.time(startInd)),datevec(datetime(2021,6,1,20,19,15))))<0.05
            powerSpec=10*log10(powerShifted);
            powerSpecNew=10*log10(specLin);

            plotRange=0.5;
            rangeInd=min(find((data.range./1000)>=plotRange));
            figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
            plot(powerSpecNew(rangeInd,:),'-b','LineWidth',2);
            ylim([-80 -30])

            set(gcf,'PaperPositionMode','auto')
            print([figdir,project,'_singleSpec_',datestr(data.time(startInd),'yyyymmdd_HHMMSS'),'_range_',num2str(plotRange),'km'],'-dpng','-r0');

            figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
            colormap('jet');
            hold on
            surf(1:size(powerSpec,2),data.range./1000,powerSpecNew,'EdgeColor','none');
            view(2)
            caxis([-60 -35]);
            colorbar
            plot([1,40],[plotRange,plotRange],'-r','LineWidth',2);
            ylabel('Range (km)');
            ylim([0 ylimUpper]);
            xlim([1 size(powerSpec,2)])
            title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))

            if flipYes
                set(gca, 'YDir','reverse');
            end

            set(gcf,'PaperPositionMode','auto')
            print([figdir,project,'_waterfall_',datestr(data.time(startInd),'yyyymmdd_HHMMSS')],'-dpng','-r0');
            stopHere=1;
        end

        %% Next beam
        startInd=endInd+1;
        ii=ii+1;
    end

    vel=vel.*data.lambda/(4*pi*prt);

    %% Plot

    disp('Plotting ...');

    f1 = figure('Position',[200 500 1800 1300],'DefaultAxesFontSize',12,'visible',showPlot);

    colormap jet

    s1=subplot(4,2,1);

    hold on
    surf(timeBeams,data.range./1000,powerV,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-110 -40]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Power (dB)')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s2=subplot(4,2,2);

    hold on
    surf(timeBeams,data.range./1000,vel,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-5 5]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Velocity (m s^{-1})')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s3=subplot(4,2,3);

    hold on
    surf(timeBeams,data.range./1000,dbz,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-60 20]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s4=subplot(4,2,4);

    hold on
    surf(timeBeams,data.range./1000,width,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([0 2]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Spectrum width (m s^{-1})')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s5=subplot(4,2,5);

    hold on
    surf(timeBeams,data.range./1000,snr,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-20 70]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Signal to noise ratio (dB)')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s6=subplot(4,2,6);

    hold on
    surf(timeBeams,data.range./1000,skew,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-1 1]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Skew (dB)')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    s8=subplot(4,2,8);

    hold on
    surf(timeBeams,data.range./1000,kurt,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([0 30]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Kurtosis (dB)')

    if flipYes
        set(gca, 'YDir','reverse');
    end

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_momentsSpec_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end