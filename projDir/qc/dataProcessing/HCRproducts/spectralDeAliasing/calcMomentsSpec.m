function moments = calcMomentsSpec(specDB,phaseVec,rx_gain,prt,lambda,noiseLev,range,dbz1km)
% Calculate moments

moments=[];

specLin=10.^(specDB./10);

sumSpecLin=sum(specLin,2,'omitnan');
sumSpecPhase=sum(specLin.*phaseVec,2,'omitnan');
sumSpecPhase2=sum(specLin.*phaseVec.^2,2,'omitnan');

% DBM
powerLin=mean(specLin,2);
moments.powerDB=10*log10(powerLin)-rx_gain;

% VEL
meanK=sumSpecPhase./sumSpecLin;
moments.vel=lambda/(4*pi*prt)*meanK;

% WIDTH
varK=(sumSpecPhase2./sumSpecLin)-meanK.^2;
sdevK=sqrt(varK);
sdevK(varK<=0)=0.0001;
moments.width=sdevK;

% SNR
noiseLin=10.^(noiseLev./10);
snrLin=(powerLin-noiseLin)./noiseLin;
snrLin(snrLin<0)=nan;
moments.snr=10*log10(snrLin);

% DBZ
range(range<0)=nan;
moments.dbz=moments.snr+20*log10(range./1000)+dbz1km;

end