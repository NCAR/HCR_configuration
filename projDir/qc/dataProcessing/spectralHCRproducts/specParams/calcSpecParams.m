function momentsSpecSmoothCorrOne=calcSpecParams(powerS,specVel,peaks1,peaks2,noiseFloor,ii,momentsSpecSmoothCorrOne)

% Find the major peak if there are two
maxPeak=peaks1;
peakDiff=peaks2(:,2)-peaks1(:,2);
maxPeak(peakDiff>0,:)=peaks2(peakDiff>0,:);

peak1Inds=find(any(~isnan(peaks1),2));
peak2Inds=find(any(~isnan(peaks2),2));

powerSS=powerS(peak1Inds,:);
specVelS=specVel(peak1Inds,:);
peaks1S=peaks1(peak1Inds,:);
peaks2S=peaks2(peak1Inds,:);
noiseFloorS=noiseFloor(peak1Inds,:);
maxPeakS=maxPeak(peak1Inds,:);

peak2IndsSS=find(any(~isnan(peaks2S),2));

% Find left edge index
[~,lindC]=max(~isnan(powerSS),[],2);

% Find right edge index
[~,rindC]=max(~isnan(fliplr(powerSS)),[],2);
rindC=size(powerSS,2)-rindC+1;

indR=1:size(powerSS,1);

lind=sub2ind(size(powerSS),indR',lindC);
rind=sub2ind(size(powerSS),indR',rindC);

maxPeaksInd=sub2ind(size(powerSS),indR',maxPeakS(:,1));
peaks1sInd=sub2ind(size(powerSS),indR',peaks1S(:,1));
peaks2sInd=sub2ind(size(powerSS),indR(peak2IndsSS)',peaks2S(peak2IndsSS,1));

% Left to right edge width
momentsSpecSmoothCorrOne.lrwidth(peak1Inds,ii)=specVelS(rind)-specVelS(lind);
% Left and right slope
momentsSpecSmoothCorrOne.lslope(peak1Inds,ii)=(noiseFloorS-maxPeakS(:,2))./(specVelS(lind)-specVelS(maxPeaksInd));
momentsSpecSmoothCorrOne.rslope(peak1Inds,ii)=(noiseFloorS-maxPeakS(:,2))./(specVelS(rind)-specVelS(maxPeaksInd));
% Left and right edge velocity
momentsSpecSmoothCorrOne.level(peak1Inds,ii)=specVelS(lind);
momentsSpecSmoothCorrOne.revel(peak1Inds,ii)=specVelS(rind);
% Left and right peak velocity
momentsSpecSmoothCorrOne.lpvel(peak1Inds,ii)=specVelS(peaks1sInd);
momentsSpecSmoothCorrOne.rpvel(peak2Inds,ii)=specVelS(peaks2sInd);
end