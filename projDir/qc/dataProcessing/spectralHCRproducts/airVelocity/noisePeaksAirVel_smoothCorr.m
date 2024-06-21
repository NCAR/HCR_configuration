function [powerRMnoiseAvRM,powerRMnoise,powerRMnoiseAvRMS,velOut,velOutS]=noisePeaksAirVel_smoothCorr(specDB,velIn,data,widthC,figdir,plotTime)

powerRMnoiseAvRM=nan(size(specDB));
powerRMnoise=nan(size(specDB));
powerRMnoiseAvRMS=nan(size(specDB));
velOut=nan(size(specDB));
%velOutS=nan(size(specDB));


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
halfSN=round(sampleNum/2);

duplicateSpec=9;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

largeInds=1:length(velSpecLarge);

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
filterAt=8;
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

    % Fill output variables
    powerRMnoiseAvRM(ii,:)=newSpecLarge(getIndsStart:getIndsStart+sampleNum-1);

    powerOrig=powerSpecLarge(ii,getIndsStart:getIndsStart+sampleNum-1);
    powOrigRMnoise=powerOrig;
    powOrigRMnoise(isnan(powerRMnoiseAvRM(ii,:)))=nan;
    powerRMnoise(ii,:)=powOrigRMnoise;

    powerSmooth=filteredSpecLarge(getIndsStart:getIndsStart+sampleNum-1);
    powerSmooth(isnan(powerRMnoiseAvRM(ii,:)))=nan;
    powerRMnoiseAvRMS(ii,:)=powerSmooth;

    velOut(ii,:)=velSpecLarge(getIndsStart:getIndsStart+sampleNum-1);

    if ismember(ii,plotRangeInds) & ~isempty(plotTime)

        close all

        %close all
        f1=figure('Position',[200 500 1000 500],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);
        t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

        s1=nexttile(1);

        hold on
        plot(velOut(ii,:),powerOrig,'-b','LineWidth',0.5);
        l1=plot(velOut(ii,:),powerRMnoise(ii,:),'-b','LineWidth',1);
        plot(velOut(ii,:),filteredSpecLarge(getIndsStart:getIndsStart+sampleNum-1),'-g','LineWidth',1);
        l2=plot(velOut(ii,:),powerRMnoiseAvRMS(ii,:),'-g','LineWidth',2);
        plot(velOut(ii,:),wcSpecLarge(getIndsStart:getIndsStart+sampleNum-1),'-r','LineWidth',1);
        l3=plot(velOut(ii,:),powerRMnoiseAvRM(ii,:),'-r','LineWidth',2);
        l4=plot(velOut(ii,:),repmat(noiseFloorAll(ii),size(velOut(ii,:))),'-c','LineWidth',1.5);
        l5=scatter(velOut(ii,peaksOut(:,1)),peaksOut(:,2),70,'m','filled','MarkerEdgeColor','black');
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
velOutS=velOut;
end