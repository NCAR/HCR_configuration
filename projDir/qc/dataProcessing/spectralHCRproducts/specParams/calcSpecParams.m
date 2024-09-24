function momentsSpecSmoothCorrOne=calcSpecParams(powerS,specVel,peaks1,peaks2,noiseFloor,ii,momentsSpecSmoothCorrOne)
loopInds=find(any(~isnan(powerS),2));

% Find the major peak if there are two
maxPeak=peaks1;
peakDiff=peaks2(:,2)-peaks1(:,2);
maxPeak(peakDiff>0,:)=peaks2(peakDiff>0,:);

for aa=1:size(loopInds,1)
    jj=loopInds(aa); % ii is the range index

    powThis=powerS(jj,:);
    velThis=specVel(jj,:);

    % Find right edge index
    rind=find(~isnan(powThis),1,'last');
    if isempty(rind)
        rind=length(powThis);
    end

    % Left to right edge width
    momentsSpecSmoothCorrOne.lrwidth(jj,ii)=velThis(rind)-velThis(1);
    % Left and right slope
    if ~isnan(maxPeak(jj,1))
        momentsSpecSmoothCorrOne.lslope(jj,ii)=(noiseFloor(jj)-maxPeak(jj,2))/(velThis(1)-velThis(maxPeak(jj,1)));
        momentsSpecSmoothCorrOne.rslope(jj,ii)=(noiseFloor(jj)-maxPeak(jj,2))/(velThis(rind)-velThis(maxPeak(jj,1)));
    end
    % Left and right edge velocity
    momentsSpecSmoothCorrOne.level(jj,ii)=velThis(1);
    momentsSpecSmoothCorrOne.revel(jj,ii)=velThis(rind);
    % Left and right peak velocity
    if ~isnan(peaks1(jj,1))
        momentsSpecSmoothCorrOne.lpvel(jj,ii)=velThis(peaks1(jj,1));
    end
    if ~isnan(peaks2(jj,1))
        momentsSpecSmoothCorrOne.rpvel(jj,ii)=velThis(peaks2(jj,1));
    end
end
end