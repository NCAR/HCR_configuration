% Analyze HCR time series
clear all;
close all;

dataDirTS='/scr/sleet3/rsfdata/projects/spicule/hcr/time_series_netcdf/';

plotTime=datetime(2021,6,1,20,19,15); % Time at which spectra should be plotted
plotRange=0.5; % Range at which spectra should be plotted

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=7.5;

startTime=datetime(2021,6,1,20,18,0);
endTime=datetime(2021,6,1,20,20,0);

%% Load data TS
disp("Getting time series data ...");

fileListTS=makeFileList(dataDirTS,startTime+seconds(1),endTime-seconds(1),'20YYMMDDxhhmmss',1);

data=[];
data.IVc=[];
data.QVc=[];

data=readHCRts(fileListTS,data,startTime,endTime);

%% Calculate moments

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

    %% Plot
    if abs(etime(datevec(data.time(startInd)),datevec(plotTime)))<0.05
        powerSpec=10*log10(powerShifted);
        powerSpecNew=10*log10(specLin);

        plotRange=0.5;
        rangeInd=min(find((data.range./1000)>=plotRange));
        figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
        plot(powerSpecNew(rangeInd,:),'-b','LineWidth',2);
        ylim([-80 -30])

        figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
        colormap('jet');
        hold on
        surf(1:size(powerSpec,2),data.range./1000,powerSpecNew,'EdgeColor','none');
        view(2)
        caxis([-60 -35]);
        colorbar
        plot([1,40],[plotRange,plotRange],'-r','LineWidth',2);
        ylabel('Range (km)');
        xlim([1 size(powerSpec,2)])
        title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))
    end

    %% Next beam
    startInd=endInd+1;
    ii=ii+1;
end
