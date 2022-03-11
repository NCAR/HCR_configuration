function [powerAdj,specVelAdj,maxIndsTest]=specDeAlias(powerSpec,duplicateSpec,sampleNum,rangeIn,plotTimeInd,maxIndsPrev)

%% Filter

maskSpec=filterGates(powerSpec);
maskSpec(1:16)=0;
powerSpecFilt=powerSpec;
powerSpecFilt(maskSpec==0,:)=nan;

plotYes=0;
if plotYes
    f1 = figure('Position',[100 500 1000 1100],'DefaultAxesFontSize',12);
    colormap jet

    subplot(1,2,1)
    surf(1:size(powerSpec,2),rangeIn./1000,powerSpec,'edgecolor','none')
    view(2)
    xlim([1,size(powerSpec,2)]);

    ylim([0 10])

    xlabel('Spectrum bin')
    ylabel('Range (km)')

    caxis([-80 -25])
    colorbar

    subplot(1,2,2)
    surf(1:size(powerSpec,2),rangeIn./1000,powerSpecFilt,'edgecolor','none')
    view(2)
    xlim([1,size(powerSpec,2)]);

    ylim([0 10])

    xlabel('Spectrum bin')
    ylabel('Range (km)')

    caxis([-80 -25])
    colorbar
end

% Add spectra side by side
powerSpecLarge=repmat(powerSpecFilt,1,duplicateSpec+2);
velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1);

powerSpecSmooth=movmedian(powerSpecLarge,round(size(powerSpec,2))/5,2);
powerSpecSmooth=powerSpecSmooth(:,sampleNum+1:end-sampleNum);
powerSpecLarge=powerSpecLarge(:,sampleNum+1:end-sampleNum);

%% Find maximum

%minInds=findMinInds(powerSpecSmooth,maskSpec);
maxInds=findMaxInds(powerSpecSmooth,maskSpec);

plotYes=1;
if plotYes
    f1 = figure('Position',[100 500 1000 1100],'DefaultAxesFontSize',12);
    colormap jet
    hold on

    surf(1:length(velSpecLarge),rangeIn./1000,powerSpecSmooth,'edgecolor','none')
    view(2)
    xlim([1,length(velSpecLarge)]);

    ylim([0 10])

    xlabel('Spectrum bin')
    ylabel('Range (km)')

    caxis([-80 -25])
    colorbar

    scatter(maxInds(maskSpec==1),rangeIn(maskSpec==1)./1000,'black','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
end


%% Compare with previous

maxIndsOrig=maxInds;

maxIndsTest=maxInds;
noiseGateInds=find(maskSpec==0);
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
specVelAdj=nan(size(powerSpec));

for kk=1:size(powerSpec,1)
    try
        powerAdj(kk,:)=powerSpecLarge(kk,maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
        specVelAdj(kk,:)=velSpecLarge(maxInds(kk)-leftAdd:maxInds(kk)+rightAdd);
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