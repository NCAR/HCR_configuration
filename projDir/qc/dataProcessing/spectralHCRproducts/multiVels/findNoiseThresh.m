function [noiseThresh,meanNoise,R2]=findNoiseThresh(powLinIn,avNum,noiseThresh)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2

powLin=powLinIn;

R2=0;
stepDown=0.0000001;

while R2<1
    powLin(powLin>noiseThresh)=[];
    % xCalc=1:length(powLin);
    sampleNum=length(powLin);
    % sig2r=(sum(xCalc.^2.*powLin,'omitmissing')./sum(powLin,'omitmissing')) ...
    %     -(sum(xCalc.*powLin,'omitmissing')./sum(powLin,'omitmissing')).^2;
    % sigN2=sampleNum.^2/12;
    meanNoise=sum(powLin,'omitmissing')/sampleNum;
    Q=sum(powLin.^2/sampleNum,'omitmissing')-meanNoise.^2;
    % R1=sigN2/sig2r;
    R2=meanNoise.^2/(Q*avNum);
    noiseThresh=noiseThresh-stepDown;
    %plot(powLin)
end
noiseThresh=10*log10(noiseThresh+stepDown);
meanNoise=10*log10(meanNoise);
end