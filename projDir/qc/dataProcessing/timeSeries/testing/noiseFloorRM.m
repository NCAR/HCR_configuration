function [sigWidthCorrRMnoise,noiseThresh,peaksOut]=noiseFloorRM(sigWidthCorr,noiseStd,sigPeaks,sigValleys,testVel,fakeMeanVel)
peakInds=find(sigPeaks==1);
if sigWidthCorr(2)<sigWidthCorr(1) & sigWidthCorr(end-1)<sigWidthCorr(end)
    if sigWidthCorr(1)>sigWidthCorr(end)
        peakInds=cat(2,peakInds,1);
    else
        peakInds=cat(2,peakInds,length(sigWidthCorr));
    end
end
peakVals=sigWidthCorr(peakInds);
pIV=cat(2,peakInds',peakVals');

peaksOut=pIV;

valleyInds=find(sigValleys==1);
valleyVals=sigWidthCorr(valleyInds);
vIV=cat(2,valleyInds',valleyVals');

vIVs=sortrows(vIV,2,'descend');

testValley=min(sigWidthCorr)+noiseStd;

if ~isempty(vIVs)
    valDiff=1;
    ii=1;
    while valDiff>0 & ii<=size(vIVs,1)
        valDiff=vIVs(ii,2)-testValley;
        if valDiff>0
            pIV(pIV(:,2)>vIVs(ii,2),:)=nan;
            vIVs(ii,:)=nan;
        end
        ii=ii+1;
    end

    vIVs(any(isnan(vIVs),2),:)=[];
    pIV(any(isnan(pIV),2),:)=[];
end

maxInd=find(pIV(:,2)==max(sigWidthCorr));

pIV(maxInd,:)=[];

% if isempty(pIV)
%     pIV=[1,sigWidthCorr(1),;length(sigWidthCorr),sigWidthCorr(end)];
% end

noiseThresh=max(pIV(:,2));

if isempty(noiseThresh)
    noiseThresh=testValley;
end

sigWidthCorrRMnoise=sigWidthCorr;
sigWidthCorrRMnoise(sigWidthCorrRMnoise<=noiseThresh)=nan;

% Check for two data stretches
rmNoiseInds=find(~isnan(sigWidthCorrRMnoise));
testDiff=diff(rmNoiseInds);
if max(abs(testDiff))>1
    [~,absMaxInd]=max(sigWidthCorr);
    if isnan(sigWidthCorrRMnoise(1)) | isnan(sigWidthCorrRMnoise(end))
        sigMask=~isnan(sigWidthCorrRMnoise);
        regs=bwconncomp(sigMask);
        for jj=1:regs.NumObjects
            thisPix=regs.PixelIdxList{jj};
            if ~ismember(absMaxInd,thisPix)
                sigWidthCorrRMnoise(thisPix)=nan;
            end
        end
    end
end

%% Return peaks
rmNoiseInds2=find(~isnan(sigWidthCorrRMnoise));
peaksOutInds=ismember(peaksOut(:,1),rmNoiseInds2);
peaksOut(peaksOutInds==0,:)=[];

end