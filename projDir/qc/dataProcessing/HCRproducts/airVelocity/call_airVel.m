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

showPlot='on';
ylimUpper=12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir,'figsAirVel/'];

file=[dataDir,infile(1:8),'/',infile];

%% Radar variables

freq=9.440617e+10;
c=299792458;
lambda=c/freq;
prt=0.000101376;
nyq=7.8311;
dbz1km_v=-23.8657;

%% Read data

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
%velFilteredTDall=nan(size(data.range,1),beamNum);
velDeAliasedTDall=nan(size(data.range,1),beamNum);
velDeAliasedSDall=nan(size(data.range,1),beamNum);
traceReflAll=nan(size(data.range,1),beamNum);
velAirAll=nan(size(data.range,1),beamNum);

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

    %% Calculate vel in time domain

    R1=mean(cIQv(:,1:end-1).*conj(cIQv(:,2:end)),2);
    velTD=lambda/(4*pi*prt)*angle(R1);

    velTDall(:,ii)=velTD;

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
    end

    % Add to output
    velDeAliasedTDall(:,ii)=velDeAliased;

    %% FFT and spectra
    
    fftIQ=fft(cIQv,[],2);

    powerRealIn=real(fftIQ).^2;
    powerImagIn=imag(fftIQ).^2;
    powerSignal=powerRealIn+powerImagIn;

    powerShifted=fftshift(powerSignal,2);

    % If nadir, reverse to get positive down
    if data.elevation(startInd)<0
        powerShifted=fliplr(powerShifted);
    end

    powerSpec=10*log10(powerShifted);

    %% De-alias in spectral domain

    if min(velMasked,[],'omitnan')<-5
        stp=1;
    end

    %prtThis=mean(prt(startInd:endInd));
    prtThis=prt;

    [powerAdj,phaseAdj]=specPowerDeAlias(powerSpec,velDeAliased,sampleNum,prtThis,lambda,data.range,velMasked);

    %% De-aliased velocity in spectral domain

    specLin=10.^(powerAdj./10);

    sumSpecLin=sum(specLin,2,'omitnan');
    sumSpecPhase=sum(specLin.*phaseAdj,2,'omitnan');

    meanK=sumSpecPhase./sumSpecLin;
    velSpecDeAliased=lambda/(4*pi*prt).*meanK;

    velDeAliasedSDall(:,ii)=velSpecDeAliased;

    %% Air velocity

    [velAir,traceRefl]=getAirVel(powerAdj,phaseAdj,sampleNum,lambda,prtThis,data.range,dbz1km_v);

    traceReflAll(:,ii)=traceRefl;
    velAirAll(:,ii)=velAir;

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

%% Plot air vel
close all

f1 = figure('Position',[200 500 1000 1200],'DefaultAxesFontSize',12,'visible',showPlot);
colM=colormap(velCols);
colormap(colM);

s1=subplot(4,1,1);
surf(timeBeams,asl./1000,velTDall,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-16 16]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
box on
title('Velocity time domain (m s^{-1})')

s2=subplot(4,1,2);
surf(timeBeams,asl./1000,velDeAliasedTDall,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-16 16]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
box on
title('Velocity time domain de-aliased (m s^{-1})')

s3=subplot(4,1,3);
surf(timeBeams,asl./1000,velDeAliasedSDall,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-16 16]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
box on
title('Velocity spectral domain de-aliased (m s^{-1})')

s4=subplot(4,1,4);
surf(timeBeams,asl./1000,velAirAll,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-16 16]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
box on
title('Air velocity de-aliased (m s^{-1})')

linkaxes([s1 s2 s3 s4],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_airVel_',datestr(timeBeams(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeBeams(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
