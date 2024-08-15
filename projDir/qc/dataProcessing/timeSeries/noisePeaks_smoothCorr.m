function [powerOrig,powerOrigRMnoise,powerSmooth,powerSmoothCorr,velOut,noiseFloorAll]=noisePeaks_smoothCorr(specDB,velIn,data,widthC,aircVel,figdir,plotTime)
powerOrig=nan(size(specDB));
powerOrigRMnoise=nan(size(specDB));
powerSmooth=nan(size(specDB));
powerSmoothCorr=nan(size(specDB));

powerSmoothAll=nan(size(specDB));
powerSmoothCorrAll=nan(size(specDB));
velOut=nan(size(specDB));

% Decide if and what to plot
plotAll=0; % Set to 1 if everything should be plotted. Plots won't be saved.
showPlot='off';

if plotAll
    plotRangeInds=18:10:size(specDB,1);
    plotTime=1;
else
    plotRangeInds=20:20:size(specDB,1);
end

sampleNum=length(data.time);

duplicateSpec=9;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

noiseFloorAll=nan(size(specDB,1),1);

%% Remove noise
% Moving average
meanOverPoints=round(sampleNum/2);
movAv=movmedian(powerSpecLarge,meanOverPoints,2);

movAv(:,1:meanOverPoints)=nan;
movAv(:,end-meanOverPoints:end)=nan;

[~,minIndTest]=min(movAv,[],2);

loopInds=find(any(~isnan(specDB),2));

testPow=nan(size(specDB,1),sampleNum);
testVel=nan(size(specDB,1),sampleNum);

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index
    testPow(ii,:)=powerSpecLarge(ii,minIndTest(ii):minIndTest(ii)+sampleNum-1);
    testVel(ii,:)=velSpecLarge(minIndTest(ii):minIndTest(ii)+sampleNum-1);
end

% Mean velocity
noiseLinV=10.^(data.noise_v./10);
sigInLin=10.^(testPow./10)-noiseLinV;

% VEL
fakeMeanVel=sum(sigInLin.*testVel,2,'omitmissing')./sum(sigInLin,2,'omitmissing');

% Filter and correct for aircraft width
filterAt=round(0.00022396*aircVel^2-0.10542*aircVel+18.132);
[sigWidthCorr,sigFiltered]=smoothAircraftWidthCorr(filterAt,testPow,fakeMeanVel,widthC,testVel,sampleNum);

% Peaks and valleys
sigPeaks=islocalmax(sigWidthCorr,2);
sigValleys=islocalmin(sigWidthCorr,2);

% Find noise floor and put spectra back together
noiseStd=5.52;

minIndTest(minIndTest>sampleNum)=minIndTest(minIndTest>sampleNum)-sampleNum;

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    % Find noise floor and remove noise
    [sigWidthCorrRMnoise,noiseFloorAll(ii),peaksOut]=noiseFloorRM(sigWidthCorr(ii,:),noiseStd,sigPeaks(ii,:),sigValleys(ii,:),testVel(ii,:),fakeMeanVel(ii));
    %[sigWidthCorrRMnoise,noiseFloorAll(ii),peaksOut]=noiseFloorRM_HilSek(sigWidthCorr(ii,:),noiseStd,sigPeaks(ii,:),sigValleys(ii,:),testVel(ii,:),fakeMeanVel(ii),aa);

    if all(isnan(sigWidthCorrRMnoise))
        continue
    end

    % Create new large spectrum
    newSpecLarge=repmat(sigWidthCorrRMnoise,1,duplicateSpec);
    
    % Check for two data stretches
    rmNoiseInds=find(~isnan(sigWidthCorrRMnoise));
    testDiff=diff(rmNoiseInds);
    if max(abs(testDiff))>1
        if ~isnan(sigWidthCorrRMnoise(1)) & ~isnan(sigWidthCorrRMnoise(end))
            firstNan=find(isnan(sigWidthCorrRMnoise),1,'first');
            newSpecLarge(1:firstNan)=nan;
        else
            stopHere=1;
        end
    end
    newSpecLarge=cat(2,nan(1,minIndTest(ii)),newSpecLarge);
    newSpecLarge(end-minIndTest(ii)+1:end)=[];

    firstPowInd=find(~isnan(newSpecLarge),1,'first');

    peaksOut(:,1)=peaksOut(:,1)-firstPowInd+minIndTest(ii)+1;
    peaksOut(peaksOut(:,1)<0,1)=peaksOut(peaksOut(:,1)<0,1)+sampleNum;

    filteredSpecLarge=repmat(sigFiltered(ii,:),1,duplicateSpec);
    filteredSpecLarge=cat(2,nan(1,minIndTest(ii)),filteredSpecLarge);
    filteredSpecLarge(end-minIndTest(ii)+1:end)=[];

    wcSpecLarge=repmat(sigWidthCorr(ii,:),1,duplicateSpec);
    wcSpecLarge=cat(2,nan(1,minIndTest(ii)),wcSpecLarge);
    wcSpecLarge(end-minIndTest(ii)+1:end)=[];

    newSpecLargeCut=newSpecLarge;
    newSpecLargeCut(1:firstPowInd-1)=[];

    specVelLargeCut=velSpecLarge;
    specVelLargeCut(1:firstPowInd-1)=[];

    % Find indices of maxima
    maxSpecLarge=max(newSpecLargeCut,[],2,'omitmissing');
    maxInds=find(newSpecLargeCut==maxSpecLarge);

    % Find correct velocity
    velAtInds=specVelLargeCut(maxInds);
    [~,velDiffMin]=min(abs(velAtInds-velIn(ii)));
    velDiffMinInd=maxInds(velDiffMin);

    sampleMult=floor(velDiffMinInd/sampleNum);
    getIndsStart=sampleMult*sampleNum+firstPowInd;

    velOut(ii,:)=velSpecLarge(getIndsStart:getIndsStart+sampleNum-1);

    powerSmoothCorr(ii,:)=newSpecLarge(getIndsStart:getIndsStart+sampleNum-1);

    % Velocity test
    y=10.^(powerSmoothCorr(ii,:)./10);
    velSpec=sum(y.*velOut(ii,:),2,'omitnan')./sum(y,2,'omitnan');

    if abs(velIn(ii)-velSpec)>1
        stop1=1;
    %     powerSmoothCorr(ii,:)=nan;
    %     velOut(ii,:)=nan;
    %     continue
    end

    % Fill output variables
    powerSmoothAll(ii,:)=filteredSpecLarge(getIndsStart:getIndsStart+sampleNum-1);
    powerSmoothCorrAll(ii,:)=wcSpecLarge(getIndsStart:getIndsStart+sampleNum-1);
    
    powerOrig(ii,:)=powerSpecLarge(ii,getIndsStart:getIndsStart+sampleNum-1);
    powOrigRMnoiseOne=powerOrig(ii,:);
    powOrigRMnoiseOne(isnan(powerSmoothCorr(ii,:)))=nan;
    powerOrigRMnoise(ii,:)=powOrigRMnoiseOne;

    powerSmoothRMnoise=powerSmoothAll(ii,:);
    powerSmoothRMnoise(isnan(powerSmoothCorr(ii,:)))=nan;
    powerSmooth(ii,:)=powerSmoothRMnoise;

    if ismember(ii,plotRangeInds) & ~isempty(plotTime)
        close all

        %close all
        f1=figure('Position',[200 500 1000 500],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

        s1=nexttile(1);

        hold on
        plot(velOut(ii,:),powerOrig(ii,:),'-b','LineWidth',0.5);
        l1=plot(velOut(ii,:),powerOrigRMnoise(ii,:),'-b','LineWidth',1);
        plot(velOut(ii,:),powerSmoothAll(ii,:),'-g','LineWidth',1);
        l2=plot(velOut(ii,:),powerSmooth(ii,:),'-g','LineWidth',2);
        plot(velOut(ii,:),powerSmoothCorrAll(ii,:),'-r','LineWidth',1);
        l3=plot(velOut(ii,:),powerSmoothCorr(ii,:),'-r','LineWidth',2);
        l4=plot(velOut(ii,:),repmat(noiseFloorAll(ii),size(velOut(ii,:))),'-c','LineWidth',1.5);
        l5=scatter(velOut(ii,peaksOut(:,1)),peaksOut(:,2),70,'m','filled','MarkerEdgeColor','black');
        ylims=s1.YLim;
        plot([velIn(ii),velIn(ii)],ylims,'-k','LineWidth',2);
        plot([velSpec,velSpec],ylims,'-m','LineWidth',2);
        hold off

        xlim([velOut(ii,1),velOut(ii,end)]);

        legend([l1,l2,l3,l4,l5],{'Original','Filtered','Filtered width corr','Noise floor','Peaks'}, ...
            'Location','northoutside','Orientation','horizontal');

        grid on
        box on

        if ~plotAll
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,'spectra/spectra_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
        end
    end
end
end