% Analyze HCR time series
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

fileList={'20210621_015305_-89.93_353.61.nc';
    '20210529_191100_-89.99_229.66.nc';
    '20210601_194538_-89.95_85.15.nc';
    '20210621_015305_-89.93_353.61.nc';
    '20210621_015437_-89.78_307.29.nc';
    '20210621_015840_89.94_315.84.nc'};

flipYes=1;

showPlot='on';

outFreq=10; % Desired output frequency in Hz

timeSpan=1/outFreq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsTS/'];

%% Loop through files
for jj=4:length(fileList)
    infile=fileList{jj};

    disp(infile);
    file=[dataDir,infile(1:8),'/',infile];

    %% Read data

    data=[];
    data.IVc=[];
    %data.IHc=[];
    data.QVc=[];
    %data.QHc=[];

    data=readHCRts(data,file);

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

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd));%./sqrt(sampleNum);

        %prtThis=mean(prt(startInd:endInd));
        prtThis=mode(data.prt);

        cIQ=cIQv;%.*sqrt(size(cIQv,2));

        R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
        R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);
        R2=mean(cIQ(:,1:end-2).*conj(cIQ(:,3:end)),2);
        R3=mean(cIQ(:,1:end-3).*conj(cIQ(:,4:end)),2);
        R4=mean(cIQ(:,1:end-4).*conj(cIQ(:,5:end)),2);
        
        powerV(:,ii)=10*log10(R0)-data.rx_gain_v;
        vel(:,ii)=data.lambda/(4*pi*prtThis)*angle(R1);
        width(:,ii)=data.lambda/(2*pi*prtThis*6^.5)*abs(log(abs(R1./R2))).^0.5;
        skew(:,ii)=abs(log(abs(R3./(R2.^3))));
        kurt(:,ii)=abs(log(abs(R4./(R2.^2))));

        % SNR
        noiseLin=10.^(data.noise_v./10);
        snrLin=(R0-noiseLin)./noiseLin;
        snrLin(snrLin<0)=nan;
        snr(:,ii)=10*log10(snrLin);
        
        % DBZ
        data.range(data.range<0)=nan;
        dbz(:,ii)=snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;
        
        timeBeams=[timeBeams;data.time(startInd)];

        startInd=endInd+1;
        ii=ii+1;
    end
   
    %% Plot

    disp('Plotting ...');

    ylimUpper=5.2;

    close all

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
    caxis([0 35]);
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
    caxis([0 20]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Kurtosis (dB)')
    
    if flipYes
        set(gca, 'YDir','reverse');
    end

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_moments_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

end