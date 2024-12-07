% Correct spectrum for broadening due to aircraft motion
function [sigWidthCorr,sigFiltered]=smoothAircraftWidthCorr(filterAt,signalIn,meanVel,corrFactor,xVel,sampleNum)

% Gaussian fit of correction signal
yWC=exp(-0.5.*((xVel-meanVel)./corrFactor).^2);

%% IFFT

ifftY=ifft(signalIn,[],2);

ifftYC=ifft(yWC,[],2);
ifftYC=ifftYC+10^(-7);

%% Correction

yC=ifftY./abs(ifftYC);
for ii=1:length(filterAt)
    yC(ii,filterAt(ii):end-filterAt(ii)+2)=0;

    %% Original uncorrected but filtered

    ifftY(ii,filterAt(ii):end-filterAt(ii)+2)=0;
end

%% FFT back
fftY=fft(yC,[],2);
sigWidthCorr=real(fftY);
sigWidthCorr=sigWidthCorr./sqrt(sampleNum).*corrFactor.*sqrt(2.*pi./((xVel(:,end)-xVel(:,1)).*(xVel(:,2)-xVel(:,1))));

sigFiltered=real(fft(ifftY,[],2));

sigWidthCorr(corrFactor<0.1,:)=sigFiltered(corrFactor<0.1,:);

end