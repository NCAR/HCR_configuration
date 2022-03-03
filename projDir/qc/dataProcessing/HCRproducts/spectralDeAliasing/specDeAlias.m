function [powerAdj,phaseAdj,maxIndsTest]=specDeAlias(powerSpec,duplicateSpec,sampleNum,rangeIn,plotTimeInd,maxIndsPrev)

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

noiseGates=sum(powerSpecFilt,2,'omitnan');
noiseMask=noiseGates~=0;
noiseMask=bwareaopen(noiseMask,2);
noiseGateInds=find(noiseMask==0);

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

%% Compare with previous

maxIndsMed=movmedian(maxIndsPrev,50,'omitnan');
maxIndsPrev=fillmissing(maxIndsPrev,'nearest');
maxIndsPrev(isnan(maxIndsMed))=nan;

maxIndsOrig=maxInds;

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

diffMaxInds=maxIndsTest-maxIndsPrev;

% Unfold
maxFolding=floor(max(abs(diffMaxInds))/sampleNum);

for ii=1:maxFolding
    maxInds(diffMaxInds>(ii-1)*sampleNum+sampleNum/2)=maxIndsOrig(diffMaxInds>(ii-1)*sampleNum+sampleNum/2)-ii*sampleNum;
    maxInds(diffMaxInds<-(ii-1)*sampleNum-sampleNum/2)=maxIndsOrig(diffMaxInds<-(ii-1)*sampleNum-sampleNum/2)+ii*sampleNum;
end

plotYes=0;
if plotYes
    plot(maxIndsTest)
    hold on
    plot(maxIndsPrev)
end

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

% Find outliers
medMaxInds=movmedian(maxIndsTest,100,'omitnan'); % 15

diffMed=maxIndsTest-medMaxInds;

maxFolding2=floor(max(abs(diffMed))/sampleNum);

for ii=1:maxFolding2
    maxInds(diffMed>(ii-1)*sampleNum+sampleNum/2)=maxIndsOrig(diffMed>(ii-1)*sampleNum+sampleNum/2)-ii*sampleNum;
    maxInds(diffMed<-(ii-1)*sampleNum-sampleNum/2)=maxIndsOrig(diffMed<-(ii-1)*sampleNum-sampleNum/2)+ii*sampleNum;
end

maxInds(diffMed>sampleNum/2 & isnan(maxIndsPrev))=maxInds(diffMed>sampleNum/2 & isnan(maxIndsPrev))-sampleNum;
maxInds(diffMed<-sampleNum/2 & isnan(maxIndsPrev))=maxInds(diffMed<-sampleNum/2 & isnan(maxIndsPrev))+sampleNum;
% maxInds(diffMed>sampleNum/2)=maxInds(diffMed>sampleNum/2)-sampleNum;
% maxInds(diffMed<-sampleNum/2)=maxInds(diffMed<-sampleNum/2)+sampleNum;

maxIndsMask=ones(size(maxInds));

maxInds(abs(diffMed)>sampleNum/4)=nan;
maxIndsMask(abs(diffMed)>sampleNum/4)=0;

maxIndsFill=maxInds;
maxIndsFill=fillmissing(maxIndsFill,'linear');

maxInds(maxIndsMask==0)=round(maxIndsFill(maxIndsMask==0));

if plotYes
    maxIndsPlot=maxInds;
    maxIndsPlot(noiseGateInds)=nan;
    scatter(1:length(maxInds),maxIndsPlot,'green')
    scatter(find(abs(diffMed)>sampleNum/4),maxIndsOrig(abs(diffMed)>sampleNum/4),'red')
    hold off
    xlim([1 400])
    ylim([2000 7000])
end

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

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

powerAdj(noiseGateInds,:)=nan;

%% Plot waterfall

if ~isempty(plotTimeInd) & ii==plotTimeInd
    close all
    %plotSpec(data,sampleNum,duplicateSpec,startInd,powerSpecLarge,ylimUpper,powerSpecFilt,showPlot,figdir,saveWaterfall)
    plotSpec(data,sampleNum,duplicateSpec,startInd,double(outRegions),ylimUpper,powerSpecFilt,showPlot,figdir,saveWaterfall)
end

end