function [sigWidthCorrRMnoise,noiseThresh]=noiseFloorRM(sigWidthCorr,noiseStd,sigPeaks,sigValleys,testVel,fakeMeanVel)
peakInds=find(sigPeaks==1);
peakVals=sigWidthCorr(peakInds);
pIV=cat(2,peakInds',peakVals');

valleyInds=find(sigValleys==1);
valleyVals=sigWidthCorr(valleyInds);
vIV=cat(2,valleyInds',valleyVals');

vIVs=sortrows(vIV,2,'descend');

testValley=min(valleyVals)+noiseStd;

valDiff=1;
ii=1;
while valDiff>0
    valDiff=vIVs(ii,2)-testValley;
    if valDiff>0
        % Find closest peaks
        valMinPeak=vIVs(ii,1)-pIV(:,1);
        prevInd=max(find(valMinPeak>0));
        nextInd=min(find(valMinPeak<0));
        vIVs(ii,:)=nan;
        pIV(prevInd,:)=nan;
        pIV(nextInd,:)=nan;
    end
    ii=ii+1;
end

vIVs(any(isnan(vIVs),2),:)=[];
pIV(any(isnan(pIV),2),:)=[];

maxInd=find(pIV(:,2)==max(sigWidthCorr));
pIV(maxInd,:)=[];

noiseThresh=max(pIV(:,2));

sigWidthCorrRMnoise=sigWidthCorr;
sigWidthCorrRMnoise(sigWidthCorrRMnoise<=noiseThresh)=nan;

% % Test velocity
% [~,velInd]=min(abs(testVel-fakeMeanVel));
% 
% rmNoiseInds=find(~isnan(sigWidthCorrRMnoise));
% 
% if ~ismember(velInd,rmNoiseInds)
%     sigWidthCorrRMnoise(:)=nan;
% end

end