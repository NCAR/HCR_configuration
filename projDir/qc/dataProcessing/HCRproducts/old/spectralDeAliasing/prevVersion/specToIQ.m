function [outputArg1,outputArg2] = specToIQ(inputArg1,inputArg2)
fftIQ=fft(cIQ,[],2)./sqrt(size(cIQ,2));
powerRealIn=real(fftIQ).^2;
powerImagIn=imag(fftIQ).^2;
powerSignal=powerRealIn+powerImagIn;
power=mean(powerSignal,2);
powerV=10*log10(power)-rx_gain;


%powerInterpDB=fftshift(powerV);
powerInterp=10.^(powerV./10);

powerRatio=powerInterp./powerOrig;
magRatio=sqrt(powerRatio);

fftRatio=fftIQ.*magRatio;

cIQfiltered(ii,:)=ifft(fftRatio).*sqrt(length(fftRatio));
end