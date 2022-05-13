% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

fileList={'20210601_174228_89.97_114.90.nc';
    '20210601_174258_68.97_346.39.nc';
    '20210601_175407_89.86_255.60.nc';
    '20210601_175438_89.95_266.91.nc';
    '20210601_175508_89.89_350.73.nc';
    '20210601_175811_89.94_214.01.nc'};

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsBacklobe/'];

%% Loop through files
for jj=4:length(fileList)
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

    %velTDall=nan(size(data.range,1),beamNum);
    dbzTDall=nan(size(data.range,1),beamNum);
    widthTDall=nan(size(data.range,1),beamNum);
    
    %% Loop through beams
    while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

        % Find start and end indices for beam
        timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
        [minDiff,endInd]=min(abs(timeDiff+timeSpan));

        sampleNum=endInd-startInd+1;

        timeBeams=[timeBeams;data.time(startInd)];
        elevBeams=[elevBeams;data.elevation(startInd)];

        prt=mode(data.prt(startInd:endInd));
        %% Window

        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

        %% Calculate vel and refle in time domain

        cIQ=cIQv.*sqrt(size(cIQv,2));

        R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
        R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);
        R2=mean(cIQ(:,1:end-2).*conj(cIQ(:,3:end)),2);

%         % VEL
%         velTD=data.lambda./(4.*pi.*prt).*angle(R1);
%         velTDall(:,ii)=velTD;

        % SNR
        noiseLin=10.^(data.noise_v./10);
        snrLin=(R0-noiseLin)./noiseLin;
        snrLin(snrLin<0)=nan;
        snr=10*log10(snrLin);
        
        % DBZ
        data.range(data.range<0)=nan;
        dbz=snr+20*log10(data.range./1000)+data.dbz1km_v;
        %dbz=10*log10(R0)-data.noise_v+data.dbz1km_v+20.*log10(data.range./1000);
        dbzTDall(:,ii)=dbz;

        % WIDTH
        width=data.lambda/(2*pi*prt*6^.5)*abs(log(abs(R1./R2))).^0.5;
        widthTDall(:,ii)=width;

        %% FFT and spectra

        fftIQ=fft(cIQv,[],2);

        powerRealIn=real(fftIQ).^2;
        powerImagIn=imag(fftIQ).^2;
        powerSignal=powerRealIn+powerImagIn;

        powerShifted=fftshift(powerSignal,2);
        
        powerSpec=10*log10(powerShifted);

        %% Plot
        close all
       
        if abs(etime(datevec(data.time(startInd)),datevec(datetime(2021,6,1,17,55,1))))<0.05

            plotRange=0.52;
            rangeInd=min(find((data.range./1000)>=plotRange));
            figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
            plot(powerSpec(rangeInd,:),'-b','LineWidth',2);
            ylim([-80 -30])

            figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
            colormap('jet');
            hold on
            surf(1:size(powerSpec,2),data.range./1000,powerSpec,'EdgeColor','none');
            view(2)
            caxis([-60 -35]);
            colorbar
            plot([1,40],[plotRange,plotRange],'-r','LineWidth',2);
            ylabel('Range (km)');
            ylim([0 ylimUpper]);
            xlim([1 size(powerSpec,2)])
            title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))
            stopHere=1;
        end
        %% Next beam
        startInd=endInd+1;
        ii=ii+1;
    end

    %% Fake asl
    asl=nan(size(dbzTDall));
    downInd=find(elevBeams<0);
    upInd=find(elevBeams>=0);
    elevIn=elevBeams';
    rangeIn=repmat(data.range,1,length(elevBeams));
    asl(:,downInd)=-1*((rangeIn(:,downInd).*cosd(abs(elevIn(downInd))-90))-10000);
    asl(:,upInd)=rangeIn(:,upInd).*cosd(abs(elevIn(upInd))-90);

    %% Plot dbz
    close all

    f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');

    s1=subplot(3,1,1);
    surf(timeBeams,asl./1000,dbzTDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([-60 0]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Reflectivity time domain (dB)')

    s2=subplot(3,1,2);
    surf(timeBeams,asl./1000,widthTDall,'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    caxis([0 3]);
    ylim([0 ylimUpper]);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    box on
    title('Width time domain (m/s)')
    
    %linkaxes([s1 s2 s3],'xy')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_dbz_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
  
end