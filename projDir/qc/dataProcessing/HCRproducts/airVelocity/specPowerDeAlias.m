function [powerAdj,specVelAdj,maxIndsTest]=specPowerDeAlias(powerSpec,velDeAlias,sampleNum,prt,lambda,rangeIn,velMasked)

%% Filter

duplicateSpec=5;

% Add spectra side by side
powerSpecLarge=repmat(powerSpec,1,duplicateSpec+2);
powerSpecSmoothLarge=movmedian(powerSpecLarge,round(size(powerSpec,2))/5,2);
powerSpecSmoothLarge=powerSpecSmoothLarge(:,sampleNum+1:end-sampleNum);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1);

powerSpecLarge=powerSpecLarge(:,sampleNum+1:end-sampleNum);
%velSpecLarge=velSpecLarge(:,sampleNum+1:end-sampleNum);

powerSpecSmooth=powerSpecSmoothLarge(:,floor(duplicateSpec/2)*sampleNum+1:end-floor(duplicateSpec/2)*sampleNum);

%% Find maximum

%minInds=findMinInds(powerSpecSmooth,maskSpec);
maxInds=findMaxIndsSpec(powerSpecSmooth,velMasked);

maxInds=maxInds+floor(duplicateSpec/2)*sampleNum;

%% Adjust maximum based on de-alias mask

phaseMax=nan(size(maxInds));

for ii=1:length(maxInds)
    if ~isnan(maxInds(ii))
        phaseMax(ii)=velSpecLarge(maxInds(ii));
    end
end
phaseDeAliased=4*pi*prt*velDeAlias/lambda;

phaseDiff=phaseDeAliased-phaseMax;

deAliasMask=zeros(size(phaseDiff));
checkFold=[2,4,6];

for jj=1:3
    deAliasMask(phaseDiff>checkFold(jj)*pi-pi)=jj;
    deAliasMask(phaseDiff<-(checkFold(jj)*pi+pi))=-jj;
end

adjMax=maxInds+deAliasMask*sampleNum;

plotYes=0;
if plotYes
    f1 = figure('Position',[100 500 1000 1100],'DefaultAxesFontSize',12);
    colormap jet
    hold on

    surf(1:size(powerSpecLarge,2),rangeIn./1000,powerSpecLarge,'edgecolor','none')
    view(2)
    xlim([1,size(powerSpecLarge,2)]);

    ylim([0 10])

    xlabel('Spectrum bin')
    ylabel('Range (km)')

    caxis([-80 -25])
    colorbar

    scatter(maxInds(~isnan(velMasked)),rangeIn(~isnan(velMasked))./1000,'black','filled');
    scatter(adjMax(~isnan(velMasked)),rangeIn(~isnan(velMasked))./1000,'red','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
end

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
        powerAdj(kk,:)=powerSpecLarge(kk,adjMax(kk)-leftAdd:adjMax(kk)+rightAdd);
        specVelAdj(kk,:)=velSpecLarge(adjMax(kk)-leftAdd:adjMax(kk)+rightAdd);
    end
end

powerAdj(isnan(velMasked),:)=nan;

%% Plot waterfall

if plotYes
    f1 = figure('Position',[100 500 500 1100],'DefaultAxesFontSize',12);
    colormap jet
    hold on

    surf(1:size(powerAdj,2),rangeIn./1000,powerAdj,'edgecolor','none')
    view(2)
    xlim([1,size(powerAdj,2)]);

    ylim([0 10])

    xlabel('Spectrum bin')
    ylabel('Range (km)')

    caxis([-80 -25])
    colorbar
end
end