function momentsSpec=calcMomentsSpec(cIQv,cIQh,sampleNum,ii,momentsSpec,data)
%% FFT and spectra

fftIQv=fft(cIQv,[],2);
fftIQh=fft(cIQh,[],2);

powerRealInV=real(fftIQv).^2;
powerImagInV=imag(fftIQv).^2;
powerRealInH=real(fftIQh).^2;
powerImagInH=imag(fftIQh).^2;

powerSignalV=powerRealInV+powerImagInV;
powerSignalH=powerRealInH+powerImagInH;

powerShiftedV=fftshift(powerSignalV,2);
powerShiftedH=fftshift(powerSignalH,2);

% Reverse to get pointing direction consistent
specLinOrigV=fliplr(powerShiftedV);
specLinOrigH=fliplr(powerShiftedH);

% Find peak
specLinV=nan(size(specLinOrigV));
specVelVecV=nan(size(specLinOrigV));
specVelVecOrig=-3*pi:2*pi/(sampleNum):3*pi;
specVelVecOrig=specVelVecOrig(1:end-1);

for kk=1:size(specLinOrigV,1)
    [~,maxInd]=max(specLinOrigV(kk,:),[],'omitnan');
    maxInd=maxInd+sampleNum;
    sBsSpec=repmat(specLinOrigV(kk,:),1,3);

    try
        specLinV(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
        specVelVecV(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2));
    catch
        specLinV(kk,:)=sBsSpec(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
        specVelVecV(kk,:)=specVelVecOrig(maxInd-floor(sampleNum/2):maxInd+floor(sampleNum/2)-1);
    end
end

%%%%%%%%%%%%%%%%%%

% DBM
powerLinV=mean(specLinOrigV,2);
momentsSpec.powerV(:,ii)=10*log10(powerLinV)-data.rx_gain_v;

powerLinH=mean(specLinOrigH,2);
momentsSpec.powerH(:,ii)=10*log10(powerLinH)-data.rx_gain_h;

% SNR
noiseLinV=10.^(data.noise_v./10);
snrLinV=(powerLinV-noiseLinV)./noiseLinV;
snrLinV(snrLinV<0)=nan;
momentsSpec.snr(:,ii)=10*log10(snrLinV);

noiseLinH=10.^(data.noise_h./10);
snrLinH=(powerLinH-noiseLinH)./noiseLinH;
snrLinH(snrLinH<0)=nan;
snrH=10*log10(snrLinH);

% DBZ
data.range(data.range<0)=nan;
momentsSpec.dbz(:,ii)=momentsSpec.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;
dbzH=snrH+20*log10(data.range./1000)+data.dbz1km_h;

% LDR
momentsSpec.ldr(:,ii)=momentsSpec.dbz(:,ii)-dbzH;

% VEL
momentsSpec.vel(:,ii)=sum(specLinV.*specVelVecV,2,'omitnan')./sum(specLinV,2,'omitnan');

% WIDTH
momentsSpec.width(:,ii)=(sum(specLinV.*(specVelVecV-momentsSpec.vel(:,ii)).^2,2,'omitnan')./sum(specLinV,2,'omitnan')).^0.5;

% SKEWNESS
momentsSpec.skew(:,ii)=sum(specLinV.*(specVelVecV-momentsSpec.vel(:,ii)).^3,2,'omitnan')./(sum(specLinV,2,'omitnan').*momentsSpec.width(:,ii).^3);

% KURTOSIS
momentsSpec.kurt(:,ii)=sum(specLinV.*(specVelVecV-momentsSpec.vel(:,ii)).^4,2,'omitnan')./(sum(specLinV,2,'omitnan').*momentsSpec.width(:,ii).^4);

end