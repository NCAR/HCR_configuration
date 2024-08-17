% Correct spectrum for broadening due to aircraft motion
function [err,errCat,sigWidthCorr,sigFiltered,signalIn1,signalIn2,sigFiltered1,sigFiltered2,inds1,inds2]=smoothingTest(filterAt,signalIn,meanVel,corrFactor,xVel,sampleNum,err)

inds1=1:2:length(signalIn);
signalIn1=signalIn(inds1);

inds2=2:2:length(signalIn);
signalIn2=signalIn(inds2);

if length(signalIn1)~=length(signalIn2)
    inds1(end)=[];
    signalIn1(end)=[];
end

% Gaussian fit of correction signal
yWC=exp(-0.5.*((xVel-meanVel)./corrFactor).^2);

%% IFFT

ifftY=ifft(signalIn,[],2);
ifftY1=ifft(signalIn1,[],2);
ifftY2=ifft(signalIn2,[],2);

ifftYC=ifft(yWC,[],2);
ifftYC=ifftYC+10^(-7);

%% Correction
notZero=2:250;

yC=ifftY./abs(ifftYC);
for ii=1:length(filterAt)
    yC(ii,filterAt(ii):end-filterAt(ii)+2)=0;
end

%% Original uncorrected but filtered

ifftYuf=repmat(ifftY,length(notZero),1,length(filterAt));
ifftYuf1=repmat(ifftY1,length(notZero),1,length(filterAt));
ifftYuf2=repmat(ifftY2,length(notZero),1,length(filterAt));
for ii=1:length(notZero)
    ifftYuf(ii,notZero(ii):end-notZero(ii)+2,:)=0;
    ifftYuf1(ii,notZero(ii):end-notZero(ii)+2,:)=0;
    ifftYuf2(ii,notZero(ii):end-notZero(ii)+2,:)=0;
end

%% FFT back
fftY=fft(yC,[],2);
sigWidthCorr=real(fftY);
sigWidthCorr=sigWidthCorr./sqrt(sampleNum).*corrFactor*sqrt(2*pi./((xVel(end)-xVel(1)).*(xVel(2)-xVel(1))));

sigFiltered=real(fft(ifftYuf,[],2));
sigFiltered1=real(fft(ifftYuf1,[],2));
sigFiltered2=real(fft(ifftYuf2,[],2));

%% Error

err12=squeeze(rmse(sigFiltered1,signalIn2,2));
err21=squeeze(rmse(sigFiltered2,signalIn1,2));

errCat=cat(2,err12,err21);

try
    err=cat(2,err,errCat);
end

end