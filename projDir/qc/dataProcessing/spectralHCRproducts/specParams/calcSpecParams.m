function momentsSpecSmoothCorrOne=calcSpecParams(powerS,specVel,peaks1,peaks2,noiseFloor,ii,momentsSpecSmoothCorrOne)
loopInds=find(any(~isnan(powerS),2));

maxPeak=peaks1;
peakDiff=peaks2(:,2)-peaks1(:,2);
maxPeak(peakDiff>0,:)=peaks2(peakDiff>0,:);

for aa=1:size(loopInds,1)
    jj=loopInds(aa); % ii is the range index

    powThis=powerS(jj,:);
    velThis=specVel(jj,:);

    rind=find(isnan(powThis),1)-1;

    momentsSpecSmoothCorrOne.lrwidth(jj,ii)=velThis(rind)-velThis(1);
    momentsSpecSmoothCorrOne.lslope(jj,ii)=(noiseFloor(jj)-maxPeak(jj,2))/(velThis(1)-velThis(maxPeak(jj,1)));
    momentsSpecSmoothCorrOne.rslope(jj,ii)=(noiseFloor(jj)-maxPeak(jj,2))/(velThis(rind)-velThis(maxPeak(jj,1)));
end
end