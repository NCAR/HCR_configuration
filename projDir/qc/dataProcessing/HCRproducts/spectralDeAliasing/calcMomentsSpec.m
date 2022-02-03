function moments = calcMomentsSpec(specDB,phaseVec,rx_gain,prt,lambda,noiseLev,range,dbz1km)
% Calculate moments

moments=[];

specLin=10.^(specDB./10);
powerLin=mean(specLin,2);

% DBM
moments.powerDB=10*log10(powerLin)-rx_gain;

% VEL
[maxPeak,maxInd]=max(specLin,[],2,'omitnan');
moments.vel=lambda/(4*pi*prt)*phaseVec(maxInd);

% WIDTH
moments.width=lambda/(2*pi*prt*6^.5)*abs(log(abs(std(specLin,0,2)))).^0.5;

%cIQ=cIQ.*sqrt(size(cIQ,2));

% %R0=mean(real(cIQ).^2+imag(cIQ).^2,2);
% R1=mean(cIQ(:,1:end-1).*conj(cIQ(:,2:end)),2);
% R2=mean(cIQ(:,1:end-2).*conj(cIQ(:,3:end)),2);
% 
% 
% % VEL
%moments.vel=lambda/(4*pi*prt)*angle(specLin);
% % WIDTH
% moments.width=lambda/(2*pi*prt*6^.5)*abs(log(abs(R1./R2))).^0.5;

% SNR
noiseLin=10.^(noiseLev./10);
snrLin=(powerLin-noiseLin)./noiseLin;
snrLin(snrLin<0)=nan;
moments.snr=10*log10(snrLin);

% DBZ
range(range<0)=nan;
moments.dbz=moments.snr+20*log10(range./1000)+dbz1km;

end