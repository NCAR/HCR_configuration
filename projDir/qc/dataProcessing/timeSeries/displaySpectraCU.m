% Analyze HCR time series
clear all;
close all;

% Data directory and file
dataDirTS='/scr/sleet3/rsfdata/projects/spicule/hcr/time_series_netcdf/20210601/';
dataFile='20210601_201905_-89.92_47.29.nc';

% Plot time and range
plotTime=datetime(2021,6,1,20,19,15); % Time at which waterfall should be plotted
plotRange=0.4; % Range at which spectra should be plotted in km

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

%% Read HCR time series data
disp('Getting time series data ...');

file=[dataDirTS,dataFile];

% Time
baseTime=ncread(file,'base_time');
timeOffset=ncread(file,'time_offset');
fileStartTime=datetime(1970,1,1)+seconds(baseTime);
data.time=fileStartTime+seconds(timeOffset)';

% Range
data.range=ncread(file,'range');

% I/Q
data.IVc=ncread([dataDirTS,dataFile],'IVc');
data.QVc=ncread([dataDirTS,dataFile],'QVc');

% PRT
data.prt=ncread(file,'prt')';
data.prt=mode(data.prt);

% Wavelength
data.lambda=ncreadatt(file,'/','radar_wavelength_cm')/100;

%% Create beams
disp('Looping through beams ...')

startInd=1;
endInd=1;
ii=1;

while endInd<=size(data.IVc,2) & startInd<size(data.IVc,2)

    % Find start and end indices for beam
    timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
    [minDiff,endInd]=min(abs(timeDiff+timeSpan));

    % Number of samples
    sampleNum=endInd-startInd+1;

    % Window data
    win=window(@hamming,sampleNum);  % Default window is Hamming
    winWeight=sampleNum/sum(win);
    winNorm=win*winWeight;

    % Create I/Q
    cIQv=winNorm'.*(data.IVc(:,startInd:endInd)+i*data.QVc(:,startInd:endInd))./sqrt(sampleNum);

    %% FFT and spectra

    fftIQ=fft(cIQv,[],2);

    powerRealIn=real(fftIQ).^2;
    powerImagIn=imag(fftIQ).^2;
    powerSignal=powerRealIn+powerImagIn;

    powerShifted=fftshift(powerSignal,2);

    % Reverse to get sign correct
    specLin=fliplr(powerShifted);

    % Calculate velocity
    specVelVec=-pi:2*pi/(sampleNum):pi;
    specVelVec=specVelVec(1:end-1);

    vel=data.lambda/(4*pi*data.prt)*specVelVec;

    %% Plot at the specified time
    if abs(etime(datevec(data.time(startInd)),datevec(plotTime)))<0.05
        % Power in dB
        powerSpecDB=10*log10(specLin);

        % Find index of specified range
        rangeInd=min(find((data.range./1000)>=plotRange));

        % Create figure
        figure('Position',[200 500 800 1200],'DefaultAxesFontSize',12);

        % Plot spectra at specified range
        s1=subplot(3,1,1);
        plot(vel,powerSpecDB(rangeInd,:),'-b','LineWidth',2);
        xlabel('Velocity (m s^{-1})');
        ylabel('Power (dB)');
        xlim([vel(1) vel(end)]);
        title({[datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS')];['Power at ',num2str(plotRange),' km range.']});
        
        % Waterfall plot
        s2=subplot(3,1,2:3);
        colormap('jet');
        hold on
        surf(vel,data.range./1000,powerSpecDB,'EdgeColor','none');
        view(2)
        caxis([-60 -30]);
        colorbar
        xlabel('Velocity (m s^{-1})');
        ylabel('Range (km)');
        xlim([vel(1) vel(end)])
        ylim([data.range(1)./1000,data.range(end)./1000])
        title('Power (dB)');

        % Plot arrow at specified range
        text(vel(1)-1.5,plotRange,'\rightarrow','FontSize',20,'color','b','FontWeight','bold','VerticalAlignment','middle');

        % Align x axis of subplots
        drawnow;
        s1Pos=s1.Position;
        s2Pos=s2.Position;
        s1.Position=[s1Pos(1),s1Pos(2),s2Pos(3),s1Pos(4)];

        % End while loop
        break
    end

    %% Next beam
    startInd=endInd+1;
    ii=ii+1;
end
