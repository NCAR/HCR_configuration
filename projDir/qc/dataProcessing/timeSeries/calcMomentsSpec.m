function momentsSpec=calcMomentsSpec(cIQv,sampleNum,ii,momentsSpec,data)
%% FFT and spectra

fftIQ=fft(cIQv,[],2);

powerRealIn=real(fftIQ).^2;
powerImagIn=imag(fftIQ).^2;
powerSignal=powerRealIn+powerImagIn;

powerShifted=fftshift(powerSignal,2);

% Reverse to get pointing direction consistent
specLinOrig=fliplr(powerShifted);

% Find peak
specLin=nan(size(specLinOrig));
specVelVec=nan(size(specLinOrig));
specVelVecOrig=-3*pi:2*pi/(sampleNum):3*pi;
specVelVecOrig=specVelVecOrig(1:end-1);

for kk=1:size(specLinOrig,1)
    [~,maxInd]=max(specLinOrig(kk,:),[],'omitnan');
    maxInd=maxInd+sampleNum;
    sBsSpec=repmat(specLinOrig(kk,:),1,3);

    try
        specLin(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
        specVelVec(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
    catch
        specLin(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
        specVelVec(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
    end
end

%%%%%%%%%%%%%%%%%%

% DBM
powerLin=mean(specLin,2);
momentsSpec.powerV(:,ii)=10*log10(powerLin)-data.rx_gain_v;

% SNR
noiseLin=10.^(data.noise_v./10);
snrLin=(powerLin-noiseLin)./noiseLin;
snrLin(snrLin<0)=nan;
momentsSpec.snr(:,ii)=10*log10(snrLin);

% DBZ
data.range(data.range<0)=nan;
momentsSpec.dbz(:,ii)=momentsSpec.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

% VEL
momentsSpec.vel(:,ii)=sum(specLin.*specVelVec,2,'omitnan')./sum(specLin,2,'omitnan');

% WIDTH
momentsSpec.width(:,ii)=(sum(specLin.*(specVelVec-momentsSpec.vel(:,ii)).^2,2,'omitnan')./sum(specLin,2,'omitnan')).^0.5;

% SKEWNESS
momentsSpec.skew(:,ii)=sum(specLin.*(specVelVec-momentsSpec.vel(:,ii)).^3,2,'omitnan')./(sum(specLin,2,'omitnan').*momentsSpec.width(:,ii).^3);

% KURTOSIS
momentsSpec.kurt(:,ii)=sum(specLin.*(specVelVec-momentsSpec.vel(:,ii)).^4,2,'omitnan')./(sum(specLin,2,'omitnan').*momentsSpec.width(:,ii).^4);

end