% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

infile='20210620_225138_89.92_169.63.nc';

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=6;
plotTimeInd=130;

plotGates=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figs/'];

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

momentsOrig.powerDB=nan(size(data.range,1),beamNum);
momentsOrig.vel=nan(size(data.range,1),beamNum);
momentsOrig.width=nan(size(data.range,1),beamNum);
momentsOrig.snr=nan(size(data.range,1),beamNum);
momentsOrig.dbz=nan(size(data.range,1),beamNum);

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
    powerSpec=10*log10(powerShifted);

    % Add spectra side by side
    powerSpecLarge=cat(2,powerSpec,powerSpec,powerSpec,powerSpec,powerSpec);

    %% Filter

    [powerSpecFilt,powerSpecMed,powerSpecMed2]=filterPowerSpec(powerSpecLarge,sampleNum);

    %% Plot waterfall

    if ii==plotTimeInd
        plotSpec(data,sampleNum,startInd,powerSpecLarge,ylimUpper,powerSpecFilt,powerSpecMed,powerSpecMed2,plotGates,showPlot,figdir)
    end

    %% Moments
    %prtThis=mean(prt(startInd:endInd));
    prtThis=prt;
    
    momentsO=calcMoments(cIQv,rx_gain_v,prtThis,lambda,noise_v,data.range,dbz1km_v);
       
    momentsOrig.powerDB(:,ii)=momentsO.powerDB;
    momentsOrig.vel(:,ii)=momentsO.vel;
    momentsOrig.width(:,ii)=momentsO.width;
    momentsOrig.snr(:,ii)=momentsO.snr;
    momentsOrig.dbz(:,ii)=momentsO.dbz;

    timeBeams=[timeBeams;data.time(startInd)];

    startInd=endInd+1;
    ii=ii+1;
end


%% Plot moments

disp('Plotting moments ...');

plotMoments('momentsOrig',momentsOrig,showPlot,timeBeams,data.range,ylimUpper,figdir,project);
