function [powerAdj,phaseAdj]=specDeAlias(powerSpec,duplicateSpec,sampleNum,rangeIn,plotTimeInd)

% Add spectra side by side
powerSpecLarge=repmat(powerSpec,1,duplicateSpec);
phaseVecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
phaseVecLarge=phaseVecLarge(1:end-1);

%% Filter

[powerSpecFilt]=filterPowerSpecPerc(powerSpecLarge,sampleNum);
powerSpecFilt(rangeIn<0,:)=nan;

%% Find spectrum boundaries

maskF=~isnan(powerSpecFilt);
maskF=bwareaopen(maskF,500,4);

powerSpecFilt(maskF==0)=nan;

distV=double(bwdist(maskF));

distV(:,1:round(sampleNum/4))=nan;
distV(:,end-round(sampleNum/4):end)=nan;

distFilt=filterDistPerc(distV,sampleNum);
distBW=~isnan(distFilt);

[leftInds,rightInds,outRegions]=getSpecBounds(distBW,powerSpecFilt,sampleNum,duplicateSpec);

%% Find maximum index between left and right
maxInds=nan(size(leftInds));

for kk=1:size(powerSpec,1)
    try
        linePiece=powerSpecLarge(kk,:);
        linePiece(1:leftInds(kk)-1)=nan;
        linePiece(rightInds(kk)+1:end)=nan;
        [maxNotUsed,maxPiece]=max(linePiece,[],'omitnan');
        maxInds(kk)=maxPiece;
    end
end

%% Build adjusted spectrum

leftAdd=floor(sampleNum/2);
if mod(sampleNum,2)==0
    rightAdd=leftAdd-1;
else
    rightAdd=leftAdd;
end

powerAdj=nan(size(powerSpec));
phaseAdj=nan(size(powerSpec));

for kk=1:size(powerSpec,1)
    try
        powerAdj(kk,:)=powerSpecLarge(kk,maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
        phaseAdj(kk,:)=phaseVecLarge(maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
    end
end

%% Plot waterfall

if ~isempty(plotTimeInd) & ii==plotTimeInd
    close all
    %plotSpec(data,sampleNum,duplicateSpec,startInd,powerSpecLarge,ylimUpper,powerSpecFilt,showPlot,figdir,saveWaterfall)
    plotSpec(data,sampleNum,duplicateSpec,startInd,double(outRegions),ylimUpper,powerSpecFilt,showPlot,figdir,saveWaterfall)
end

end