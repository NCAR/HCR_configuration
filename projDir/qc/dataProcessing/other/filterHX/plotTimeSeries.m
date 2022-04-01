% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec

% Input file
%infile='20210625_193014_-90.00_220.66.nc';
%infile='20210625_193518_-89.97_114.69.nc';
infile='20210625_194027_-89.98_105.11.nc';

% Time for waterfall plot
plotTime=datetime(2021,6,25,19,40,30);
% Range for spectral plot
plotRange=5.35; % km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

%figdir=[dataDir,'testHX/'];

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

file=[dataDir,infile(1:8),'/',infile];

disp(infile);

data=[];
data.IHc=[];
data.QHc=[];

data=readHCRts(data,file);

%% Prepare processing
beamNum=ceil(size(data.IHc,2)/(timeSpan*10000));

powerDB=nan(size(data.range,1),beamNum);

timeBeams=[];
elevBeams=[];

startInd=1;
endInd=1;
ii=1;

% Loop through beams
while endInd<=size(data.IHc,2) & startInd<size(data.IHc,2)

    % Find start and end indices for beam
    timeDiff=etime(datevec(data.time(startInd)),datevec(data.time));
    [minDiff,endInd]=min(abs(timeDiff+timeSpan));

    sampleNum=endInd-startInd+1;

    timeBeams=[timeBeams;data.time(startInd)];
    elevBeams=[elevBeams;data.elevation(startInd)];

    % Window
    win=window(@hamming,sampleNum);  % Default window is Hamming
    winWeight=sampleNum/sum(win);
    winNorm=win*winWeight;

    cIQh=winNorm'.*(data.IHc(:,startInd:endInd)+i*data.QHc(:,startInd:endInd))./sqrt(sampleNum);

    %% FFT and spectra

    fftIQ=fft(cIQh,[],2);

    powerRealIn=real(fftIQ).^2;
    powerImagIn=imag(fftIQ).^2;
    powerSignal=powerRealIn+powerImagIn;

    powerShifted=fftshift(powerSignal,2);
    powerSpec=10*log10(powerShifted);

    %% Plot waterfall
    if abs(etime(datevec(data.time(startInd)),datevec(plotTime)))<0.03
        f1 = figure('Position',[100 500 600 1100],'DefaultAxesFontSize',12);

        subplot(3,1,1)
        xSpec=-pi:2*pi/(sampleNum-1):pi;

        [~,rangeGateInd]=min(abs(data.range-plotRange*1000));
        plot(xSpec,powerSpec(rangeGateInd,:),'-b','LineWidth',1);
        title(['Range: ',num2str(data.range(rangeGateInd)./1000,3),' km']);

        xlim([-pi pi])
        
        colormap jet

        subplot(3,1,2:3)
        surf(xSpec,data.range./1000,powerSpec,'edgecolor','none')
        view(2)
        xlim([-pi pi])

        ylim([-0.5 15])

        xlabel('Spectrum bin')
        ylabel('Range (km)')

        title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))
        hold on
        plot([-pi,-2.8],[data.range(rangeGateInd)./1000,data.range(rangeGateInd)./1000],'-k','LineWidth',2)

        caxis([-80 -25])
        colorbar
    end

    %% Calculate power

    prtThis=prt;

    moments=[];

    cIQ=cIQh.*sqrt(size(cIQh,2));

    R0=mean(real(cIQ).^2+imag(cIQ).^2,2);

    powerDB(:,ii)=10*log10(R0)-rx_gain_h;

    startInd=endInd+1;
    ii=ii+1;
end

%% Plot power

f1=figure('Position',[700 900 800 700],'DefaultAxesFontSize',12);
colormap('jet');
surf(timeBeams,data.range/1000,powerDB,'edgecolor','none');
view(2);
ylabel('km');
caxis([-105 -92]);
ylim([-0.5 15]);
xlim([timeBeams(1),timeBeams(end)]);

hold on
plot([plotTime-seconds(0.5),plotTime+seconds(0.5)],[data.range(rangeGateInd)./1000,data.range(rangeGateInd)./1000],'-k','LineWidth',1)
plot([plotTime,plotTime],[-0.5,15],'-k','LineWidth',1)
colorbar
grid on
title('H power (dB)')

%print(f1,[figdir,project,'_H_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
