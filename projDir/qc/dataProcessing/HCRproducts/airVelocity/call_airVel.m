% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

infile='20210529_191100_-89.99_229.66.nc';
%infile='20210620_225107_83.48_16.92.nc';
%infile='20210620_225138_89.92_169.63.nc';
%infile='20210621_015305_-89.93_353.61.nc';
%infile='20210621_015437_-89.78_307.29.nc';
%infile='20210621_015840_89.94_315.84.nc';

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;
duplicateSpec=7; % Number of duplicates of spectra

showPlot='on';
ylimUpper=6;
plotTimeInd=[];
saveWaterfall=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsAirVel/'];

file=[dataDir,infile(1:8),'/',infile];

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
momentsOrigSpec.airVel=nan(size(data.range,1),beamNum);

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

    %% FFT and spectra

    fftIQ=fft(cIQv,[],2);

    powerRealIn=real(fftIQ).^2;
    powerImagIn=imag(fftIQ).^2;
    powerSignal=powerRealIn+powerImagIn;

    powerShifted=fftshift(powerSignal,2);
    powerShifted=fliplr(powerShifted);
    powerSpec=10*log10(powerShifted);

    %% De-alias

    [powerAdj,phaseAdj]=specDeAlias(powerSpec,duplicateSpec,sampleNum,data.range,plotTimeInd);

    %% Air velocity

    %prtThis=mean(prt(startInd:endInd));
    prtThis=prt;

    [momentsOrigSpec.airVel(:,ii),momentsOrigSpec.traceRefl(:,ii)]=getAirVel(powerAdj,phaseAdj,sampleNum,lambda,prtThis,data.range,dbz1km_v);

    %% Moments
    
    if ii==100
        stopHere=1;
    end

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

    timeBeams=[timeBeams;data.time(startInd)];

    startInd=endInd+1;
    ii=ii+1;
end


%% Plot moments

disp('Plotting moments ...');

plotMoments('momentsOrigIQ',momentsOrigIQ,showPlot,timeBeams,data.range,ylimUpper,figdir,project);

plotMoments('momentsOrigSpec',momentsOrigSpec,showPlot,timeBeams,data.range,ylimUpper,figdir,project);

%% Plot air vel
f1 = figure('Position',[200 500 600 1000],'DefaultAxesFontSize',12,'visible',showPlot);
colormap jet

subplot(2,1,1)
surf(timeBeams,data.range./1000,momentsOrigSpec.airVel,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-16 16]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Air velocity (m s^{-1})')

subplot(2,1,2)
surf(timeBeams,data.range./1000,momentsOrigSpec.traceRefl,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-40 30]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Tracer reflectivity (dBZ)')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_airVel_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');