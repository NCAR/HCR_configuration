% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

fileList={'20210529_191100_-89.99_229.66.nc';
    '20210601_200353_-89.98_196.46.nc';
    '20210601_201905_-89.92_47.29.nc';
    '20210621_015305_-89.93_353.61.nc';
    '20210621_015437_-89.78_307.29.nc'};

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=11;
ylimLower=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsAirVel/'];

%% Radar variables
freq=9.440617e+10;
c=299792458;
lambda=c/freq;
prt=0.000101376;
nyq=7.8311;
dbz1km_v=-23.8657;
noise_v=-61.301;

%% Loop through files
for jj=1:length(fileList)
    infile=fileList{jj};

    disp(infile);
    file=[dataDir,infile(1:8),'/',infile];
    
    %% Read data

    data=[];
    data.IVc=[];
    data.QVc=[];

    data=readHCRts(data,file);

    %% Prepare processing
    beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

    timeBeams=[];
    elevBeams=[];

    startInd=1;
    endInd=1;
    ii=1;

    velTDall=nan(size(data.range,1),beamNum);
    dbzTDall=nan(size(data.range,1),beamNum);
    %snrTDall=nan(size(data.range,1),beamNum);
    %velFilteredTDall=nan(size(data.range,1),beamNum);
    velDeAliasedTDall=nan(size(data.range,1),beamNum);
    velDeAliasedSDall=nan(size(data.range,1),beamNum);
    traceReflAll=nan(size(data.range,1),beamNum);
    velAirAll=nan(size(data.range,1),beamNum);
    velSmallerAll=nan(size(data.range,1),beamNum);
    velLargerAll=nan(size(data.range,1),beamNum);
    velMaxAll=nan(size(data.range,1),beamNum);

    %% Set up for de-aliasing

    defaultPrev=nyq;

    velPrev=repmat(defaultPrev,length(data.range),1);
    prevCount=zeros(size(velPrev));
    prevKeep=nan(size(velPrev));
    flipYes=0;

    %% Loop through beams
    while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

        % Find start and end indices for beam
        timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
        [minDiff,endInd]=min(abs(timeDiff+timeSpan));

        sampleNum=endInd-startInd+1;

        timeBeams=[timeBeams;data.time(startInd)];
        elevBeams=[elevBeams;data.elevation(startInd)];

        %% Window

        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

        %% Calculate vel and refl in time domain

        cIQ=cIQv.*sqrt(size(cIQv,2));

        R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
        R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);

        % VEL
        velTD=lambda/(4*pi*prt)*angle(R1);
        velTDall(:,ii)=velTD;

        % SNR
        noiseLin=10.^(noise_v./10);
        snrLin=(R0-noiseLin)./noiseLin;
        snrLin(snrLin<0)=nan;
        snr=10*log10(snrLin);
        %snrTDall(:,ii)=snr;

        % DBZ
        data.range(data.range<0)=nan;
        dbz=snr+20*log10(data.range./1000)+dbz1km_v;
        dbzTDall(:,ii)=dbz;

        %% Filter gates in time domain

        velMasked=filterGatesTD(velTD);
        %velFilteredTDall(:,ii)=velMasked;

        %% De-alias
        if data.elevation(startInd)>0
            velRay=-velMasked;
        else
            velRay=velMasked;
        end

        velDeAliased=deAliasSingleRay(velRay,velPrev,nyq,[],[]);

        % Set up time consistency check

        [velPrev,prevCount,prevKeep,flipYes]=setUpPrev(velDeAliased,velPrev,prevCount,prevKeep,flipYes,data.elevation(startInd),outFreq,defaultPrev);

        if data.elevation(startInd)>0
            velDeAliased=-velDeAliased;
            velRay=-velRay;
        end

        % Add to output
        velDeAliasedTDall(:,ii)=velDeAliased;

        %% FFT and spectra

        fftIQ=fft(cIQv,[],2);

        powerRealIn=real(fftIQ).^2;
        powerImagIn=imag(fftIQ).^2;
        powerSignal=powerRealIn+powerImagIn;

        powerShifted=fftshift(powerSignal,2);

        % Reverse to get pointing direction consistent
        powerShifted=fliplr(powerShifted);
        
        powerSpec=10*log10(powerShifted);

        %% De-alias in spectral domain

        if min(velMasked,[],'omitnan')<-5
            stp=1;
        end

        %prtThis=mean(prt(startInd:endInd));
        prtThis=prt;

        [powerAdj,phaseAdj]=specPowerDeAlias(powerSpec,velDeAliased,sampleNum,prtThis,lambda,data.range,velRay);

        %% De-aliased velocity in spectral domain

        specLin=10.^(powerAdj./10);

        sumSpecLin=sum(specLin,2,'omitnan');
        sumSpecPhase=sum(specLin.*phaseAdj,2,'omitnan');

        meanK=sumSpecPhase./sumSpecLin;
        velSpecDeAliased=lambda/(4*pi*prt).*meanK;

        velDeAliasedSDall(:,ii)=velSpecDeAliased;

        %% Air velocity

        [velAir,velSmaller,velLarger,velMax,traceRefl]=getAirVel(powerAdj,phaseAdj,data.elevation(startInd),sampleNum,lambda,prtThis,data.range,dbz1km_v,noise_v);

        traceReflAll(:,ii)=traceRefl;
        velAirAll(:,ii)=velAir;
        velSmallerAll(:,ii)=velSmaller;
        velLargerAll(:,ii)=velLarger;
        velMaxAll(:,ii)=velMax;

        startInd=endInd+1;
        ii=ii+1;
    end

    %% Fake asl
    asl=nan(size(velTDall));
    downInd=find(elevBeams<0);
    upInd=find(elevBeams>=0);
    elevIn=elevBeams';
    rangeIn=repmat(data.range,1,length(elevBeams));
    asl(:,downInd)=-1*((rangeIn(:,downInd).*cosd(abs(elevIn(downInd))-90))-10000);
    asl(:,upInd)=rangeIn(:,upInd).*cosd(abs(elevIn(upInd))-90);

    %% Plot vel
    close all

    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colM=colormap(velCols);
    colormap(colM);

    s1=subplot(3,1,1);
    surf(timeBeams,asl./1000,velTDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity time domain (m s^{-1})')

    s2=subplot(3,1,2);
    surf(timeBeams,asl./1000,velDeAliasedTDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity time domain de-aliased (m s^{-1})')

    s3=subplot(3,1,3);
    surf(timeBeams,asl./1000,velDeAliasedSDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity spectral domain de-aliased (m s^{-1})')

    linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_vel_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

    %% Plot air vel and reflectivities

    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');

    s1=subplot(3,1,1);
    surf(timeBeams,asl./1000,dbzTDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-40 20]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Reflectivity time domain (dBZ)');

    s2=subplot(3,1,2);
    surf(timeBeams,asl./1000,traceReflAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-40 20]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Tracer reflectivity time domain (dBZ)');

    s3=subplot(3,1,3);
    surf(timeBeams,asl./1000,velAirAll,'edgecolor','none');
    s3.Colormap=colM;
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Air velocity de-aliased (m s^{-1})')

    linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_reflAirVel_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

    %% Plot max, smaller, and larger vel
    
    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap(colM);

    s1=subplot(3,1,1);
    surf(timeBeams,asl./1000,velMaxAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity max peak (m s^{-1})')

    s2=subplot(3,1,2);
    surf(timeBeams,asl./1000,velLargerAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity larger vel peak (m s^{-1})')

    s3=subplot(3,1,3);
    surf(timeBeams,asl./1000,velSmallerAll,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-16 16]);
    ylim([ylimLower ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Velocity smaller vel peak (m s^{-1})')

    linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velMaxLargerSmaller_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end