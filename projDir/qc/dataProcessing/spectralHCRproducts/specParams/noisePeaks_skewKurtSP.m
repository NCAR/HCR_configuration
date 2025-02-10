function [powerOrig,powerOrigRMnoise,powerSmooth,powerSmoothCorr,velOut,noiseFloorAllMov,peakIndsAll1,peakIndsAll2]= ...
    noisePeaks_skewKurtSP(specDB,data,widthC,filterAt,sampleTime,figdir,plotTime)

% Initialize output
powerOrig=nan(size(specDB));
powerOrigRMnoise=nan(size(specDB));
powerSmooth=nan(size(specDB));

% Decide if and what to plot
plotAll=0; % Set to 1 if everything should be plotted. Plots won't be saved.
showPlot='on';

if plotAll
    plotRangeInds=18:1:size(specDB,1);
    plotTime=1;
else
    plotRangeInds=20:20:size(specDB,1);
end

% Create velocity vector
sampleNum=length(data.time);
velSpec=-pi:2*pi/(sampleNum):pi;
velSpec=velSpec(1:end-1).*data.lambda./(4*pi.*data.prt);

%% Noise floor
% Noise floor averaging number
avNum=3;

% Find noise floor
rawPowBig=repmat(specDB,1,3);
rawPowMov=movmean(rawPowBig,avNum,2);
rawPowMov=rawPowMov(:,sampleNum+1:2*sampleNum);

noiseFloorAll=findNoiseThreshMat(rawPowMov,avNum);

% Average noise floor
noiseFloorAllF=fillmissing(noiseFloorAll,'linear');
noiseFloorAllMov=movmedian(noiseFloorAllF,5);
noiseFloorAllMov(isnan(noiseFloorAll))=nan;

%% Aircraft width correction
% Mean velocity
sigInLin=10.^(specDB./10);

% VEL
meanVel=sum(sigInLin.*velSpec,2,'omitmissing')./sum(sigInLin,2,'omitmissing');

% Filter and correct for aircraft width
[sigWidthCorr,sigFiltered]=smoothAircraftWidthCorr(filterAt,specDB,meanVel,widthC,velSpec,sampleNum);

%% Remove spectral noise
powerSmoothCorr=sigWidthCorr;
powerSmoothCorr(powerSmoothCorr<noiseFloorAllMov)=nan;

%% Reorganize and clean
% Reorganize so full spectra fit within bounds (no folding) and remove data
% above noise that are too insignificant
sigDupl=repmat(powerSmoothCorr,1,2);
velDupl=-pi:2*pi/(sampleNum):3*pi;
velDupl=velDupl(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,2));
velOut=repmat(velSpec,size(sigDupl,1),1);

% If we are in debug mode, fill out other ouput and plot variables
if ~isempty(plotTime)
    % Fill output variables
    powerOrig=specDB;
    pOdupl=repmat(powerOrig,1,2);

    powerOrigRMnoise=powerOrig;
    pORMnoiseDupl=repmat(powerOrigRMnoise,1,2);

    powerSmooth=sigFiltered;
    pSdupl=repmat(powerSmooth,1,2);

    powerSmoothAll=sigFiltered;
    pSduplAll=repmat(powerSmoothAll,1,2);

    powerSmoothCorrAll=sigWidthCorr;
    pSCduplAll=repmat(powerSmoothCorrAll,1,2);
end

% Remove two stretches that are not valid
loopInds=find(any(~isnan(sigDupl),2));
for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    sigThis=sigDupl(ii,:);

    sigMask=isnan(sigThis);
    if sum(sigMask)==0
        [~,firstEmpty]=min(sigThis);
    else
        firstEmpty=find(isnan(sigThis),1,'first');
    end

    if firstEmpty~=1 & sum(sigMask)~=0
        sigThisDiff=diff(sigMask);
        startNan=find(sigThisDiff==1);
        endNan=find(sigThisDiff==-1);
        startNan=startNan(1:length(endNan));
        nanDiffMax=max(endNan-startNan);
        sigMaskFilled=bwareaopen(sigMask,nanDiffMax-1);
        firstEmpty=find(sigMaskFilled==1,1,'first');
    end
    if firstEmpty~=1
        powerSmoothCorr(ii,:)=sigThis(firstEmpty:firstEmpty+sampleNum-1);
        velOut(ii,:)=velDupl(firstEmpty:firstEmpty+sampleNum-1);
        if ~isempty(plotTime)
            powerOrig(ii,:)=pOdupl(ii,firstEmpty:firstEmpty+sampleNum-1);
            powerOrigRMnoise(ii,:)=pORMnoiseDupl(ii,firstEmpty:firstEmpty+sampleNum-1);
            powerSmooth(ii,:)=pSdupl(ii,firstEmpty:firstEmpty+sampleNum-1);
            powerSmoothAll(ii,:)=pSduplAll(ii,firstEmpty:firstEmpty+sampleNum-1);
            powerSmoothCorrAll(ii,:)=pSCduplAll(ii,firstEmpty:firstEmpty+sampleNum-1);
        end
    end

    % Check for two data stretches
    sigWidthCorrRMnoise=powerSmoothCorr(ii,:);
    rmNoiseInds=find(~isnan(sigWidthCorrRMnoise));    
    testDiff=diff(rmNoiseInds);
    if max(abs(testDiff))>1
        spread=max(sigWidthCorrRMnoise)-noiseFloorAllMov(ii);
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
        powerSmoothCorr(ii,:)=sigWidthCorrRMnoise;
    end
end

if ~isempty(plotTime)
    powerOrigRMnoise(isnan(powerSmoothCorr))=nan;
    powerSmooth(isnan(powerSmoothCorr))=nan;
end

%% Peaks
sigPeaks=islocalmax(powerSmoothCorr,2,'MaxNumExtrema',2);
peakIndsAll1=nan(size(powerSmoothCorr,1),2);
peakIndsAll2=nan(size(powerSmoothCorr,1),2);

noPeaks=sum(sigPeaks,2);
powerOrig(noPeaks==0,:)=nan;
powerOrigRMnoise(noPeaks==0,:)=nan;
powerSmooth(noPeaks==0,:)=nan;
powerSmoothCorr(noPeaks==0,:)=nan;
velOut(noPeaks==0,:)=nan;
noiseFloorAllMov(noPeaks==0,:)=nan;

loopInds=find(any(~isnan(powerSmoothCorr),2));
for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    peakInds=find(sigPeaks(ii,:)==1);
    peakVals=powerSmoothCorr(ii,peakInds);
    peaksOut=cat(2,peakInds',peakVals');
    
    if ~isempty(peaksOut)
        peakIndsAll1(ii,:)=peaksOut(1,:);
        if size(peaksOut,1)==2 & ~isnan(peaksOut(2,2))
            peakIndsAll2(ii,:)=peaksOut(2,:);
        end
    end

    %% Plot
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