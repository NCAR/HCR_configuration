function momentsTime=calcMomentsTime(cIQv,cIQh,ii,momentsTime,data)
cIQv=cIQv.*sqrt(size(cIQv,2));

R0v=mean(real(cIQv).^2+imag(cIQv).^2,2);
R1v=mean(cIQv(:,1:end-1).*conj(cIQv(:,2:end)),2);
R2v=mean(cIQv(:,1:end-2).*conj(cIQv(:,3:end)),2);
R3v=mean(cIQv(:,1:end-3).*conj(cIQv(:,4:end)),2);
R4v=mean(cIQv(:,1:end-4).*conj(cIQv(:,5:end)),2);

R0h=mean(real(cIQh).^2+imag(cIQh).^2,2);

momentsTime.powerV(:,ii)=10*log10(R0v)-data.rx_gain_v;
momentsTime.vel(:,ii)=data.lambda/(4*pi*mode(data.prt))*angle(R1v);
momentsTime.width(:,ii)=data.lambda/(2*pi.*mode(data.prt)*6^.5)*abs(log(abs(R1v./R2v))).^0.5;
momentsTime.skew(:,ii)=abs(log(abs(R3v./(R2v.^3))));
momentsTime.kurt(:,ii)=abs(log(abs(R4v./(R2v.^2))));

momentsTime.powerH(:,ii)=10*log10(R0h)-data.rx_gain_h;

% SNRV
noiseLinV=10.^(data.noise_v./10);
snrLinV=(R0v-noiseLinV)./noiseLinV;
snrLinV(snrLinV<0)=nan;
momentsTime.snr(:,ii)=10*log10(snrLinV);

noiseLinH=10.^(data.noise_h./10);
snrLinH=(R0h-noiseLinH)./noiseLinH;
snrLinH(snrLinH<0)=nan;
snrH=10*log10(snrLinH);

% DBZ
data.range(data.range<0)=nan;
momentsTime.dbz(:,ii)=momentsTime.snr(:,ii)+20*log10(data.range./1000)+data.dbz1km_v;

dbzH=snrH+20*log10(data.range./1000)+data.dbz1km_h;

% LDR
momentsTime.ldr(:,ii)=momentsTime.dbz(:,ii)-dbzH;

end

