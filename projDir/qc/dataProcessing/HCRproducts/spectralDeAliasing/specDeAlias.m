% Analyze HCR time series

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
freqData='dummy';
qcVersion='dummy';

infile='20210529_191100_-89.99_229.66.nc';

outFreq=10; % Desired output frequency in Hz
timeSpan=1/outFreq;

showPlot='on';
ylimUpper=6;
plotTimeInd=150;
plotRangeInd=70;

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

%% Read data

data.IVc=[];
data.IHc=[];
data.QVc=[];
data.QHc=[];

data=readHCRts(data,file);

%% Calculate moments
beamNum=ceil(size(data.IVc,2)/(timeSpan*10000));

powerV=nan(size(data.range,1),beamNum);
vel=nan(size(data.range,1),beamNum);
width=nan(size(data.range,1),beamNum);

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

    powerSpecMed=movmedian(powerSpec,size(powerSpec,2)/8,2);

    powerSpecMin=min(powerSpecMed,[],2,'omitnan');
    powerSpecFilt=powerSpec;
    powerSpecFilt(powerSpec<powerSpecMin)=nan;

    if ii==plotTimeInd
        f1 = figure('Position',[100 500 400 1100],'DefaultAxesFontSize',12,'visible',showPlot);

        colormap jet

        surf(1:sampleNum,data.range./1000,powerSpecFilt,'EdgeColor','none');
        view(2);

        xlim([1,sampleNum]);
        ylim([0,data.range(end)./1000]);

        ylim([0 ylimUpper])

        xlabel('Sample number')
        ylabel('Range (km)')
        title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))

        caxis([-100 25])
        colorbar

        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,'waterfall_',datestr(data.time(startInd),'yyyymmdd_HHMMSS')],'-dpng','-r0');
        
        f1 = figure('Position',[100 500 600 400],'DefaultAxesFontSize',12,'visible',showPlot);

        hold on
        colormap jet

        plot(powerSpec(plotRangeInd,:));
        plot(powerSpecMed(plotRangeInd,:));
        plot(powerSpecFilt(plotRangeInd,:));
        xlabel('Sample number')
        ylabel('Power (dB)')

        title([datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'),' range ',num2str(data.range(plotRangeInd)./1000,2),' km'])

        print(f1,[figdir,'spectrum_',datestr(data.time(startInd),'yyyymmdd_HHMMSS'),'_range_',num2str(data.range(plotRangeInd)./1000),'km'],'-dpng','-r0');
    end

    %prtThis=mean(prt(startInd:endInd));
    prtThis=prt;

    %    R0=mean(cIQ.^2,2);
    R1=mean(cIQv(:,1:end-1).*conj(cIQv(:,2:end)),2);
    R2=mean(cIQv(:,1:end-2).*conj(cIQv(:,3:end)),2);

    %     power=sqrt(abs(R0));
    %     powerDB(:,ii)=10*log10(power)-rx_gain;
    power=mean(powerSignal,2);
    powerV(:,ii)=10*log10(power)-rx_gain_v;

    %     R0=mean(cIQv.^2,2);
    %     R1=mean(cIQv(:,1:end-1).*conj(cIQv(:,2:end)),2);
    %     R2=mean(cIQv(:,1:end-2).*conj(cIQv(:,3:end)),2);

    %powerV(:,ii)=10*log10(sqrt(abs(R0)))-rx_gain_v;
    vel(:,ii)=lambda/(4*pi*prtThis)*angle(R1);
    width(:,ii)=lambda/(2*pi*prtThis*6^.5)*abs(log(abs(R1./R2))).^0.5;

    timeBeams=[timeBeams;data.time(startInd)];

    startInd=endInd+1;
    ii=ii+1;
end

powerV(:,all(isnan(powerV),1))=[];
vel(:,all(isnan(vel),1))=[];
width(:,all(isnan(width),1))=[];

%% Plot

disp('Plotting ...');

f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

s1=subplot(3,1,1);

hold on
surf(timeBeams,data.range./1000,powerV,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-110 -70]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Power (dB)')
s1pos=s1.Position;
s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];

s2=subplot(3,1,2);

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
s2pos=s2.Position;
s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

s3=subplot(3,1,3);

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
s3pos=s3.Position;
s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_moments_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');

