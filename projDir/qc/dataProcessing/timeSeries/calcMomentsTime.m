function momentsTime=calcMomentsTime(cIQv,prt,ii,momentsTime,data)
cIQ=cIQv.*sqrt(size(cIQv,2));

R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);
R2=mean(cIQ(:,1:end-2).*conj(cIQ(:,3:end)),2);
R3=mean(cIQ(:,1:end-3).*conj(cIQ(:,4:end)),2);
R4=mean(cIQ(:,1:end-4).*conj(cIQ(:,5:end)),2);

momentsTime.powerV(:,ii)=10*log10(R0)-data.rx_gain_v;
momentsTime.vel(:,ii)=data.lambda/(4*pi*prt)*angle(R1);
momentsTime.width(:,ii)=data.lambda/(2*pi*prt*6^.5)*abs(log(abs(R1./R2))).^0.5;
momentsTime.skew(:,ii)=abs(log(abs(R3./(R2.^3))));
momentsTime.kurt(:,ii)=abs(log(abs(R4./(R2.^2))));

% SNR
noiseLin=10.^(data.noise_v./10);
snrLin=(R0-noiseLin)./noiseLin;
snrLin(snrLin<0)=nan;
momentsTime.snr(:,ii)=10*log10(snrLin);

% DBZ
data.range(data.range<0)=nan;
momentsTime.dbz(:,ii)=momentsTime.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

end

