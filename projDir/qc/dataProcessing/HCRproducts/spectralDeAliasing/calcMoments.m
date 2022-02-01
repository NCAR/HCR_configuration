function moments = calcMoments(cIQ,rx_gain,prt,lambda,noiseLev,range,dbz1km)
% Calculate moments

moments=[];

cIQ=cIQ.*sqrt(size(cIQ,2));

R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);
R2=mean(cIQ(:,1:end-2).*conj(cIQ(:,3:end)),2);

% DBM
moments.powerDB=10*log10(R0)-rx_gain;
% VEL
moments.vel=lambda/(4*pi*prt)*angle(R1);
% WIDTH
moments.width=lambda/(2*pi*prt*6^.5)*abs(log(abs(R1./R2))).^0.5;

% SNR
noiseLin=10.^(noiseLev./10);
snrLin=(R0-noiseLin)./noiseLin;
snrLin(snrLin<0)=nan;
moments.snr=10*log10(snrLin);

% DBZ
range(range<0)=nan;
moments.dbz=moments.snr+20*log10(range./1000)+dbz1km;

end