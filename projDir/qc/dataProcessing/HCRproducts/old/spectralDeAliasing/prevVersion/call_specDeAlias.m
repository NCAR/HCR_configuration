% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

fileList={'20210529_191100_-89.99_229.66.nc';
    '20210620_225107_83.48_16.92.nc';
    '20210620_225138_89.92_169.63.nc';
    '20210620_230015_89.89_308.46.nc';
    '20210620_232102_87.04_268.26.nc';
    '20210620_233311_89.82_89.28.nc';
    '20210621_015305_-89.93_353.61.nc';
    '20210621_015437_-89.78_307.29.nc';
    '20210621_015840_89.94_315.84.nc';
    '20210621_015910_89.96_132.69.nc';
    '20210621_015941_89.89_71.42.nc';
    '20210621_035228_-89.94_332.98.nc';
    '20210624_233336_89.87_28.29.nc';
    '20210625_215734_-89.77_83.01.nc'};

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;
duplicateSpec=7; % Number of duplicates of spectra

showPlot='on';
plotTimeInd=[];
saveWaterfall=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsDeAlias/'];

%% Radar variables

freq=9.440617e+10;
c=299792458;
lambda=c/freq;
rx_gain_v=45.9;
rx_gain_h=45.5;
prt=0.000101376;
dbz1km_v=-23.8657;
noise_v=-61.301;

%% Read data

for jj=1:length(fileList)

    infile=fileList{jj};
    file=[dataDir,infile(1:8),'/',infile];

    disp(infile);

    data=[];
    data.IVc=[];
    data.IHc=[];
    data.QVc=[];
    data.QHc=[];

    data=readHCRts(data,file);

    %% Prepare processing
    beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

    momentsOrigIQ.powerDB=nan(size(data.range,1),beamNum);
    momentsOrigIQ.vel=nan(size(data.range,1),beamNum);
    momentsOrigIQ.width=nan(size(data.range,1),beamNum);
    momentsOrigIQ.snr=nan(size(data.range,1),beamNum);
    momentsOrigIQ.dbz=nan(size(data.range,1),beamNum);

    momentsOrigSpec.powerDB=nan(size(data.range,1),beamNum);
    momentsOrigSpec.vel=nan(size(data.range,1),beamNum);
    momentsOrigSpec.width=nan(size(data.range,1),beamNum);
    momentsOrigSpec.snr=nan(size(data.range,1),beamNum);
    momentsOrigSpec.dbz=nan(size(data.range,1),beamNum);

    prevIndsAll=nan(size(data.range,1),beamNum);

    timeBeams=[];
    elevBeams=[];

    startInd=1;
    endInd=1;
    ii=1;

    % Loop through beams
    while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

        % Find start and end indices for beam
        timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
        [minDiff,endInd]=min(abs(timeDiff+timeSpan));

        sampleNum=endInd-startInd+1;

        timeBeams=[timeBeams;data.time(startInd)];
        elevBeams=[elevBeams;data.elevation(startInd)];

        %disp(data.time(startInd));

        % Window
        win=window(@hamming,sampleNum);  % Default window is Hamming
        winWeight=sampleNum/sum(win);
        winNorm=win*winWeight;

        cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);
        %cIQv=(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

        %% FFT and spectra

        fftIQ=fft(cIQv,[],2);

        powerRealIn=real(fftIQ).^2;
        powerImagIn=imag(fftIQ).^2;
        powerSignal=powerRealIn+powerImagIn;

        powerShifted=fftshift(powerSignal,2);
        powerShifted=fliplr(powerShifted);
        powerSpec=10*log10(powerShifted);

        %% De-alias

        % Set up previous maximum index
        if ii==1
            prevCount=zeros(size(fftIQ,1),1);
            maxIndsPrev=nan(size(fftIQ,1),1);
            maxIndsPrev(:)=size(powerSpec,2)*duplicateSpec/2;
            maxIndsKeep=nan(size(fftIQ,1),1);
        end
        prevIndsAll(:,ii)=maxIndsPrev;

        % De-alias
        [powerAdj,phaseAdj,maxIndsPrev]=specDeAlias(powerSpec,duplicateSpec,sampleNum,data.range,plotTimeInd,maxIndsPrev);

        % Clean up max indices
        prevLarge=maxIndsPrev;
        prevLarge=movmean(prevLarge,5,'omitnan');
        prevLarge=movmean(prevLarge,5,'includenan');
        prevMask=~isnan(prevLarge);
        prevMask=bwareaopen(prevMask,10);
        maxIndsPrev(prevMask==0)=nan;
        powerAdj(prevMask==0,:)=nan;

        % Replenish old max indeces
        maxIndsKeep(~isnan(maxIndsPrev))=maxIndsPrev(~isnan(maxIndsPrev));
        prevCount(~isnan(maxIndsPrev))=0;
        prevCount(isnan(maxIndsPrev))=prevCount(isnan(maxIndsPrev))+1;

        % Add old max indeces
        maxIndsPrev(isnan(maxIndsPrev) & prevCount<outFreq*5 & ~isnan(maxIndsKeep))=maxIndsKeep(isnan(maxIndsPrev) & prevCount<outFreq*5 & ~isnan(maxIndsKeep));

        % Clean up old max indeces
        maxIndsKeep(prevCount>=outFreq*5)=nan;

        % Interpolate prev max
        maxIndsMed=movmedian(maxIndsPrev,50,'omitnan');
        maxIndsMedLarge=movmedian(maxIndsPrev,100,'omitnan');

        maxIndsPrev(isnan(maxIndsPrev) & ~isnan(maxIndsMed))=maxIndsMedLarge(isnan(maxIndsPrev) & ~isnan(maxIndsMed));

        nanIndsMask=isnan(maxIndsPrev);
        nanIndsMask=bwareaopen(nanIndsMask,50);
        nanIndsMask(~isnan(maxIndsMedLarge))=0;

        maxIndsPrev(nanIndsMask==1)=size(powerSpec,2)*duplicateSpec/2;
        maxIndsPrev=fillmissing(maxIndsPrev,'linear','EndValues','nearest');

        %% Moments
        %prtThis=mean(prt(startInd:endInd));
        prtThis=prt;

        momentsOIQ=calcMomentsIQ(cIQv,rx_gain_v,prtThis,lambda,noise_v,data.range,dbz1km_v);
       
        momentsOrigIQ.powerDB(:,ii)=momentsOIQ.powerDB;
        momentsOrigIQ.vel(:,ii)=momentsOIQ.vel;
        momentsOrigIQ.width(:,ii)=momentsOIQ.width;
        momentsOrigIQ.snr(:,ii)=momentsOIQ.snr;
        momentsOrigIQ.dbz(:,ii)=momentsOIQ.dbz;

        momentsOS=calcMomentsSpec(powerAdj,phaseAdj,rx_gain_v,prtThis,lambda,noise_v,data.range,dbz1km_v);

        momentsOrigSpec.powerDB(:,ii)=momentsOS.powerDB;
        momentsOrigSpec.vel(:,ii)=momentsOS.vel;
        momentsOrigSpec.width(:,ii)=momentsOS.width;
        momentsOrigSpec.snr(:,ii)=momentsOS.snr;
        momentsOrigSpec.dbz(:,ii)=momentsOS.dbz;

        startInd=endInd+1;
        ii=ii+1;
    end

    %% Fake asl
    asl=nan(size(momentsOrigSpec.vel));
    downInd=find(elevBeams<0);
    upInd=find(elevBeams>=0);
    elevIn=elevBeams';
    rangeIn=repmat(data.range,1,length(elevBeams));
    asl(:,downInd)=-1*((rangeIn(:,downInd).*cosd(abs(elevIn(downInd))-90))-10000);
    asl(:,upInd)=rangeIn(:,upInd).*cosd(abs(elevIn(upInd))-90);

    %% Post processing

    nyquistVel=7.8311;

    [velFinal,changeMat]=postProcessDeAlias(momentsOrigSpec.vel,nyquistVel);

    %% Get ylimits
    dbzSum=abs(sum(momentsOrigSpec.dbz,2,'omitnan'));
    maxGate=max(find(dbzSum>50));
    aslGate=median(asl(maxGate,:));
    if isempty(upInd)
        ylimits=[aslGate-500,10100];
    else
        ylimits=[0,aslGate+500];
    end
    ylimits=ylimits/1000;

    %% Plot final
    close all

    f1=figure('Position',[200 500 800 600],'DefaultAxesFontSize',12);
    colormap('jet');
    %surf(timeBeams,asl/1000,velFinal,'edgecolor','none');
    surf(timeBeams,asl/1000,momentsOrigSpec.vel,'edgecolor','none');
    view(2);
    ylabel('km');
    caxis([-15 15]);
    ylim(ylimits);
    xlim([timeBeams(1),timeBeams(end)]);
    colorbar
    grid on
    title('Velocity (m s^{-1})')
    %print(f1,[figdir,project,'_velFinal_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    print(f1,[figdir,project,'_vel_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    %% Plot moments

    disp('Plotting moments ...');

    plotMoments('momentsOrigIQ',momentsOrigIQ,showPlot,timeBeams,asl,ylimits,figdir,project);

    %plotMoments('momentsOrigSpec',momentsOrigSpec,showPlot,timeBeams,asl,ylimits,figdir,project);
end