% Correct spectrum for broadening due to aircraft motion
function [sigWidthCorr,sigFiltered]=aircraftWidthCorr(signalIn,corrFactor,xVel,noiseThresh,sampleNum)
noiseLinV=10.^(noiseThresh./10);
sigInLin=10.^(signalIn./10)-noiseLinV;

% VEL
meanVel=sum(sigInLin.*xVel,'omitmissing')/sum(sigInLin,'omitmissing');

% Gaussian fit of correction signal
yWC=exp(-0.5.*((xVel-meanVel)/corrFactor).^2);

%% IFFT

ifftY=ifft(signalIn,[],2);

ifftYC=ifft(yWC,[],2);
ifftYC=ifftYC+10^(-7);

%% Correction
notZero=6;

yC=ifftY./abs(ifftYC);
yC(notZero:end-notZero+2)=0;

%% Original uncorrected but filtered

ifftYuf=ifftY;
ifftYuf(notZero:end-notZero+2)=0;

%% FFT back
fftY=fft(yC,[],2);
sigWidthCorr=real(fftY);
sigWidthCorr=sigWidthCorr./sqrt(sampleNum).*corrFactor*sqrt(2*pi./((xVel(end)-xVel(1)).*(xVel(2)-xVel(1))));

sigFiltered=real(fft(ifftYuf,[],2));

end