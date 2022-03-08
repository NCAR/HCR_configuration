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
        [~,maxPiece]=max(linePiece,[],'omitnan');
        maxInds(kk)=maxPiece;
    end
end

%% Compare with previous

maxIndsOrig=maxInds;

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

plotYes=0;
if plotYes
    plot(maxIndsTest)
    hold on
    plot(maxIndsPrev)
    xlim([1 400])
    ylim([2000 7000])
end

diffMaxInds=maxIndsTest-maxIndsPrev;

% Unfold
maxFolding=floor(max(abs(diffMaxInds))/sampleNum);

for ii=1:maxFolding
    maxInds(diffMaxInds>(ii-1)*sampleNum+sampleNum/2)=maxIndsOrig(diffMaxInds>(ii-1)*sampleNum+sampleNum/2)-ii*sampleNum;
    maxInds(diffMaxInds<-(ii-1)*sampleNum-sampleNum/2)=maxIndsOrig(diffMaxInds<-(ii-1)*sampleNum-sampleNum/2)+ii*sampleNum;
end

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

%% Find outlier areas and correct with sample num
medMaxInds=movmedian(maxIndsTest,50,'omitnan');

diffMed=maxIndsTest-medMaxInds;

maxFolding2=floor(max(abs(diffMed))/sampleNum);

for ii=1:maxFolding2
    maxInds(diffMed>(ii-1)*sampleNum+sampleNum/2)=maxIndsOrig(diffMed>(ii-1)*sampleNum+sampleNum/2)-ii*sampleNum;
    maxInds(diffMed<-(ii-1)*sampleNum-sampleNum/2)=maxIndsOrig(diffMed<-(ii-1)*sampleNum-sampleNum/2)+ii*sampleNum;
end

%% Find individual outliers

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

medMaxInds2=movmedian(maxIndsTest,15,'omitnan');
diffMed2=maxIndsTest-medMaxInds2;

maxIndsMask=ones(size(maxInds));
maxIndsMask(abs(diffMed2)>sampleNum/4)=0;
maxIndsMask2=zeros(size(maxIndsMask));
% maxIndsMask2=~maxIndsMask;
% maxIndsMask2=bwareaopen(maxIndsMask2,5);

maxInds(maxIndsMask==0 & maxIndsMask2==0)=nan;

maxIndsTest=maxInds;
maxIndsTest(noiseGateInds)=nan;

maxIndsFill=maxIndsTest;
maxIndsFill=fillmissing(maxIndsFill,'linear');

maxInds(maxIndsMask==0 & maxIndsMask2==0)=round(maxIndsFill(maxIndsMask==0 & maxIndsMask2==0));

if plotYes
    maxIndsPlot=maxInds;
    maxIndsPlot(noiseGateInds)=nan;
    scatter(1:length(maxInds),maxIndsPlot,'green')
    scatter(find(maxIndsMask==0 & maxIndsMask2==0),maxIndsOrig(maxIndsMask==0 & maxIndsMask2==0),'red')
    hold off
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