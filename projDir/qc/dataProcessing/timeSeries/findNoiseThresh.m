function [noiseThresh,meanNoiseOut,R2]=findNoiseThresh(powIn,avNum)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2

powLin2=10.^(powIn./10);

testThresh=sort(powLin2,'descend');

testPow=testThresh;

ii=1;
R2=0;

while R2<1
    testPow(testPow>testThresh(ii))=[];
    sampleNum=length(testPow);
    meanNoise=sum(testPow)./sampleNum;
    Q=sum(testPow.^2./sampleNum)-meanNoise.^2;
    R2=meanNoise.^2./(Q*avNum);
    ii=ii+1;
end

noisePow=powLin2;
noisePow(noisePow>testThresh(ii-1))=[];
meanNoisePow=sum(noisePow)./sampleNum;
noiseThresh=real(10*log10(testThresh(ii-1)));
meanNoiseOut=10*log10(meanNoisePow);
end