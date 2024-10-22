function [powerOrig,powerOrigRMnoise,powerSmooth,powerSmoothCorr,velOut,noiseFloorAllMov,peakIndsAll1,peakIndsAll2]= ...
    noisePeaks_smoothCorr(specDB,velIn,data,widthC,aircVel,sampleTime,figdir,plotTime)
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
    plotRangeInds=18:1:size(specDB,1);
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

% Moving average
meanOverPoints=round(sampleNum/2);
movAv=movmedian(powerSpecLarge,meanOverPoints,2);

movAv(:,1:meanOverPoints)=nan;
movAv(:,end-meanOverPoints:end)=nan;

[~,minIndTest]=min(movAv,[],2);

loopInds=find(any(~isnan(specDB),2));

testPow=nan(size(specDB,1),sampleNum);
testVel=nan(size(specDB,1),sampleNum);

% % Noise floor averaging number
% avNum=3;
% noiseFloorAll=nan(size(specDB,1),1);
% 
% for aa=1:size(loopInds,1)
%     ii=loopInds(aa); % ii is the range index
%     testPow(ii,:)=powerSpecLarge(ii,minIndTest(ii):minIndTest(ii)+sampleNum-1);
%     testVel(ii,:)=velSpecLarge(minIndTest(ii):minIndTest(ii)+sampleNum-1);
% 
%     % Find noise floor
%     rawPowBig=repmat(testPow(ii,:),1,3);
%     rawPowMov=movmean(rawPowBig,avNum);
%     rawPowMov=rawPowMov(sampleNum+1:2*sampleNum);
% 
%     [noiseFloorAll(ii),~,~]=findNoiseThresh(rawPowMov,avNum);
% end

% Noise floor averaging number
avNum=3;

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index
    testPow(ii,:)=powerSpecLarge(ii,minIndTest(ii):minIndTest(ii)+sampleNum-1);
    testVel(ii,:)=velSpecLarge(minIndTest(ii):minIndTest(ii)+sampleNum-1);
end

% Find noise floor
rawPowBig=repmat(testPow,1,3);
rawPowMov=movmean(rawPowBig,avNum,2);
rawPowMov=rawPowMov(:,sampleNum+1:2*sampleNum);

noiseFloorAll=findNoiseThreshMat(rawPowMov,avNum);

%% Average noise floor
noiseFloorAllF=fillmissing(noiseFloorAll,'linear');
noiseFloorAllMov=movmedian(noiseFloorAllF,5);
noiseFloorAllMov(isnan(noiseFloorAll))=nan;

% Mean velocity
sigInLin=10.^(testPow./10);

% VEL
fakeMeanVel=sum(sigInLin.*testVel,2,'omitmissing')./sum(sigInLin,2,'omitmissing');

% Filter and correct for aircraft width
if sampleTime==0.1
    filterAt=round(0.00022396.*aircVel.^2-0.10542.*aircVel+18.132);
elseif sampleTime==0.01
    filterAt=round(0.000126.*aircVel.^2-0.0548.*aircVel+10.7);
else
    error('Sample time must be 0.1 or 0.01.')
end
filterAt=fillmissing(filterAt,'nearest');
[sigWidthCorr,sigFiltered]=smoothAircraftWidthCorr(filterAt,testPow,fakeMeanVel,widthC,testVel,sampleNum);

%Remove below noise
sigWidthCorrRMnoiseAll=sigWidthCorr;
sigWidthCorrRMnoiseAll(sigWidthCorrRMnoiseAll<noiseFloorAllMov)=nan;

% Remove two stretches that are not valid
for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    sigWidthCorrRMnoise=sigWidthCorrRMnoiseAll(ii,:);    
    if all(isnan(sigWidthCorrRMnoise))
        continue
    end

    % Check for two data stretches
    rmNoiseInds=find(~isnan(sigWidthCorrRMnoise));
    testDiff=diff(rmNoiseInds);
    if max(abs(testDiff))>1
        cutReg=0;
        spread=max(sigWidthCorrRMnoise)-noiseFloorAllMov(ii);
        if ~isnan(sigWidthCorrRMnoise(1)) & ~isnan(sigWidthCorrRMnoise(end))
            cutReg=1;
            firstNan=find(isnan(sigWidthCorrRMnoise),1);
            sigWidthCorrRMnoise=cat(2,sigWidthCorrRMnoise(firstNan:end),sigWidthCorrRMnoise(1:firstNan-1));
        end
        sigMask=~isnan(sigWidthCorrRMnoise);
        regs=bwconncomp(sigMask);
        for jj=1:regs.NumObjects
            maxReg=max(sigWidthCorrRMnoise(regs.PixelIdxList{jj}));
            spreadReg=maxReg-noiseFloorAllMov(ii);
            if sampleTime==0.1
                if spread~=spreadReg & (spreadReg<1 | spreadReg/spread<0.3)
                    sigWidthCorrRMnoise(regs.PixelIdxList{jj})=nan;
                end
            else
                if spread~=spreadReg & (spreadReg<5 | spreadReg/spread<0.3)
                    sigWidthCorrRMnoise(regs.PixelIdxList{jj})=nan;
                end
            end
        end
        if cutReg
            sigWidthCorrRMnoise=cat(2,sigWidthCorrRMnoise(end-firstNan+2:end),sigWidthCorrRMnoise(1:end-firstNan+1));
        end
        sigWidthCorrRMnoiseAll(ii,:)=sigWidthCorrRMnoise;
    end
end

% Peaks
sigPeaks=islocalmax(sigWidthCorrRMnoiseAll,2,'MaxNumExtrema',2);
peakIndsAll1=nan(size(sigWidthCorrRMnoiseAll,1),2);
peakIndsAll2=nan(size(sigWidthCorrRMnoiseAll,1),2);

minIndTest(minIndTest>sampleNum)=minIndTest(minIndTest>sampleNum)-sampleNum;
for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    sigWidthCorrRMnoise=sigWidthCorrRMnoiseAll(ii,:);    
    if all(isnan(sigWidthCorrRMnoise))
        continue
    end
    % Create new large spectrum
    newSpecLarge=repmat(sigWidthCorrRMnoise,1,duplicateSpec);

    newSpecLarge=cat(2,nan(1,minIndTest(ii)),newSpecLarge);
    newSpecLarge(end-minIndTest(ii)+1:end)=[];

    % Remove spectra pieces that are not complete
    largeMask=~isnan(newSpecLarge);
    diffLA=diff(largeMask);
    startLA=find(diffLA==1)+1;
    endLA=find(diffLA==-1);

    if ~isempty(startLA) & ~isempty(endLA)
        if endLA(1)<startLA(1)
            startLA=[1,startLA];
        end
        if endLA(end)<startLA(end)
            endLA=[endLA,length(newSpecLarge)];
        end

        areas=endLA-startLA+1;

        ua=unique(areas);
        if length(ua)>1
            countRegs=nan(1,length(ua));
            for kk=1:length(ua)
                countRegs(kk)=sum(areas==ua(kk));
            end
            remove1=find(countRegs==1);
            if ~isempty(remove1)
                for ll=1:length(remove1)
                    regInd=find(areas==ua(remove1(ll)));
                    newSpecLarge(startLA(regInd):endLA(regInd))=nan;
                end
            end
        end
    end

    firstPowInd=find(~isnan(newSpecLarge),1,'first');

    % Peaks
    peakInds=find(sigPeaks(ii,:)==1);
    if sigWidthCorr(ii,2)<sigWidthCorr(ii,1) & sigWidthCorr(ii,end-1)<sigWidthCorr(ii,end)
        if sigWidthCorr(ii,1)>sigWidthCorr(ii,end)
            peakInds=cat(2,peakInds,1);
        else
            peakInds=cat(2,peakInds,sampleNum);
        end
    end
    peakVals=sigWidthCorrRMnoise(peakInds);
    pIV=cat(2,peakInds',peakVals');
    pIV(any(isnan(pIV),2),:)=[];
    if size(pIV,1)>2
        pIV=sortrows(pIV,2,'descend');
        pIV=pIV(1:2,:);
    end

    peaksOut=pIV;
    peaksOut(:,1)=peaksOut(:,1)-firstPowInd+minIndTest(ii)+1;
    peaksOut(peaksOut(:,1)<0,1)=peaksOut(peaksOut(:,1)<0,1)+sampleNum;
    peaksOut=sortrows(peaksOut);
    if ~isempty(peaksOut)
        peakIndsAll1(ii,:)=peaksOut(1,:);
        if size(peaksOut,1)==2 & ~isnan(peaksOut(2,2))
            peakIndsAll2(ii,:)=peaksOut(2,:);
        end
    end

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
        l4=plot(velOut(ii,:),repmat(noiseFloorAllMov(ii),size(velOut(ii,:))),'-c','LineWidth',1.5);
        l5=scatter(velOut(ii,peakIndsAll1(ii,1)),peakIndsAll1(ii,2),70,'m','filled','MarkerEdgeColor','black');
        if ~isnan(peakIndsAll2(ii,1))
            scatter(velOut(ii,peakIndsAll2(ii,1)),peakIndsAll2(ii,2),70,'m','filled','MarkerEdgeColor','black');
        end
        ylims=s1.YLim;
        plot([velIn(ii),velIn(ii)],ylims,'-k','LineWidth',2);
        plot([velSpec,velSpec],ylims,'-m','LineWidth',2);
        hold off

        xlim([velOut(ii,1),velOut(ii,end)]);

        legend([l1,l2,l3,l4,l5],{'Original','Filtered','Filtered width corr','Noise floor','Peaks'}, ...
            'Location','northoutside','Orientation','horizontal');

        xlabel('Velocity (m s^{-1})');
        ylabel('Power (dB)')

        grid on
        box on

        if ~plotAll
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,'spectra/spectra_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
        end
    end
end

%% Waterfall plot

if ~isempty(plotTime) & ~all(isnan(velOut(:)))
    xlims=[floor(min(velOut(:),[],'omitmissing')),ceil(max(velOut(:),[],'omitmissing'))];

    f1=figure('Position',[200 500 1600 885],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');
    t=tiledlayout(1,4,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);
    hold on
    surf(velOut,data.range./1000,powerOrig,'edgecolor','none');
    view(2)
    xlim(xlims);
    clim([-60,-10]);
    xlabel('Velocity (m s^{-1})')
    ylabel('Range (km)')
    title('Raw')
    box on
    grid on
    ylimits=s1.YLim;
    for kk=1:length(plotRangeInds)
        plot(xlims,[data.range(plotRangeInds(kk))/1000,data.range(plotRangeInds(kk))/1000],'-m','LineWidth',1.5);
    end
    s1.YLim=ylimits;
    s1.SortMethod='childorder';

    s2=nexttile(2);
    hold on
    surf(velOut,data.range./1000,powerOrigRMnoise,'edgecolor','none');
    view(2)
    xlim(xlims);
    clim([-60,-10]);
    xlabel('Velocity (m s^{-1})')
    ylabel('Range (km)')
    title('Raw noise removed')
    box on
    grid on

    s3=nexttile(3);
    hold on
    surf(velOut,data.range./1000,powerSmooth,'edgecolor','none');
    view(2)
    xlim(xlims);
    clim([-60,-10]);
    xlabel('Velocity (m s^{-1})')
    ylabel('Range (km)')
    title('Filtered')
    box on
    grid on

    s4=nexttile(4);
    hold on
    surf(velOut,data.range./1000,powerSmoothCorr,'edgecolor','none');
    view(2)
    xlim(xlims);
    clim([-60,-10]);
    xlabel('Velocity (m s^{-1})')
    ylabel('Range (km)')
    title('Filtered and corrected')
    box on
    grid on
    for aa=1:size(loopInds,1)
        kk=loopInds(aa);
        if ~isnan(peakIndsAll1(kk,1))
            scatter(velOut(kk,peakIndsAll1(kk,1)),data.range(kk)/1000,20,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
        end
        if ~isnan(peakIndsAll2(kk,1))
            scatter(velOut(kk,peakIndsAll2(kk,1)),data.range(kk)/1000,20,'filled','MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor','k');
        end
    end
    s4.SortMethod='childorder';
    colorbar
    if ~plotAll
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,'waterfall/waterfall_',datestr(plotTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0');
    end
end
end