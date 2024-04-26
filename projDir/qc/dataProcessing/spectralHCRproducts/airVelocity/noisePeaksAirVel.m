function [powerRMnoiseAvRM,velOut,peakVelsOut,peakPowsOut]=noisePeaksAirVel(specDB,velIn,data,firFilt,filtShift,widthC,plotTime)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries

sampleNum=length(data.time);
halfSN=round(sampleNum/2);

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));
midInd=round(length(velSpecLarge)/2);
midInds=midInd-halfSN:midInd+halfSN;

largeInds=1:length(velSpecLarge);

powerRMnoiseAv=nan(size(specDB));
powerRMnoiseAvRM=nan(size(specDB));
powerRMnoiseRaw=nan(size(specDB));
powerRMnoiseRawPlot=nan(size(specDB));
velOut=nan(size(specDB));

meanNoiseAll=nan(size(specDB,1),1);
noiseThreshAll=meanNoiseAll;

%% Remove noise
% Moving average
meanOverPoints=3; % Average over this number of points
secondMean=25;
movAv=movmedian(powerSpecLarge,meanOverPoints,2);
movAv2=movmedian(movAv,secondMean,2);

movAv(:,1:round(sampleNum/10))=nan;
movAv(:,end-round(sampleNum/10):end)=nan;

movAv2(:,1:round(sampleNum/10))=nan;
movAv2(:,end-round(sampleNum/10):end)=nan;

loopInds=find(any(~isnan(specDB),2));

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    thisMov=movAv(ii,:);

    % Find noise floor and noise threshold    
    [noiseThreshAll(ii),meanNoiseAll(ii)]=findNoiseThresh(thisMov(midInds),meanOverPoints);

    % Remove noise below threshold
    thisMovRM=movAv2(ii,:);
    thisMovRM(thisMovRM<noiseThreshAll(ii))=nan;
    thisMask1=~isnan(thisMovRM);
    %thisMask=bwareafilt(thisMask1,5);
    thisMovRM(thisMask1==0)=nan;

    % Find indices of maxima
    maxMov=max(thisMov,[],2,'omitmissing');
    maxInds=find(thisMov==maxMov);

    % Find correct velocity
    velAtInds=velSpecLarge(maxInds);
    [~,velDiffMin]=min(abs(velAtInds-velIn(ii)));
    velDiffMinInd=maxInds(velDiffMin);

    smallerInds=largeInds(velDiffMinInd-halfSN:velDiffMinInd+halfSN);
    smallerMask=thisMask1(smallerInds);
    smallerMov=thisMov(smallerInds);
    maxInd1=find(smallerMov==maxMov);

    specRegs=bwconncomp(smallerMask);

    for jj=1:specRegs.NumObjects
        if ~ismember(maxInd1,specRegs.PixelIdxList{jj})
            smallerMask(specRegs.PixelIdxList{jj})=0;
        end
    end

    maskOneInd=find(smallerMask,1,'first');
    finalInds=smallerInds(maskOneInd):smallerInds(maskOneInd)+sampleNum-1;

    if ~isempty(finalInds)
        largeMask=zeros(size(largeInds));
        largeMask(smallerInds)=smallerMask;

        thisMov(largeMask==0)=nan;
        thisMov=thisMov(finalInds);

        thisMovRM(largeMask==0)=nan;
        thisMovRM=thisMovRM(finalInds);

        thisRaw=powerSpecLarge(ii,:);
        thisRawPlot=thisRaw;
        thisRaw(largeMask==0)=nan;
        thisRaw=thisRaw(finalInds);
        thisRawPlot=thisRawPlot(finalInds);

        thisVel=velSpecLarge(finalInds);

        powerRMnoiseAv(ii,:)=thisMov;
        powerRMnoiseAvRM(ii,:)=thisMovRM;
        powerRMnoiseRaw(ii,:)=thisRaw;
        powerRMnoiseRawPlot(ii,:)=thisRawPlot;
        velOut(ii,:)=thisVel;
    else
        continue
    end
end
    
%% Find peaks and plot spectra

peaksIndsAll=nan(size(velOut,1),2);
peakVelsOut=nan(size(velOut,1),2);
peakPowsOut=nan(size(velOut,1),2);

loopInds2=find(any(~isnan(powerRMnoiseAvRM),2));

% FIR
powerApp=cat(2,nan(size(powerRMnoiseAvRM,1),filtShift),powerRMnoiseAvRM,nan(size(powerRMnoiseAvRM,1),filtShift));  % Append D zeros to the input data

powFilt=filter(firFilt,powerApp');
powFilt=powFilt';
powerSmoothAll=powFilt(:,2*filtShift+1:end);

% Width Correction
for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index
    corrAW=aircraftWidthCorr(powerRMnoiseRawPlot(ii,:),widthC,velOut(ii,:),data.noise_v,noiseThreshAll(ii));
    %corrAW=aircraftWidthCorr(powerRMnoiseAv(ii,:),widthC,velOut(ii,:),data.noise_v,noiseThreshAll(ii));
    %corrAW=aircraftWidthCorr(powerRMnoiseAvRM(ii,:),widthC,velOut(ii,:),data.noise_v,noiseThreshAll(ii));
end

% Find peaks
minDiffMS=1.5; % Minimum velocity difference between peaks in m/s
minDiffPix=round(minDiffMS/(velSpecLarge(2)-velSpecLarge(1)));
findPeaks=islocalmax(powerSmoothAll,2,'MinProminence',0.5,'FlatSelection','center', ...
    'MinSeparation',minDiffPix,'MaxNumExtrema',2);

% Decide if and what to plot
plotAll=0; % Set to 1 if everything should be plotted. Plots won't be saved.
showPlot='off';

if plotAll
    plotRangeInds=18:10:size(specDB,1);
    plotTime=1;
else
    plotRangeInds=20:20:size(specDB,1);
end

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/figsMultiVel/cases/';
for aa=1:size(loopInds2,1)
    ii=loopInds2(aa); % ii is the range index

    thisPeaks=findPeaks(ii,:);
    peakInds=find(thisPeaks==1);

    peaksIndsAll(ii,1:length(peakInds))=peakInds;
    peakVelsOut(ii,1:length(peakInds))=velOut(ii,peakInds);
    peakPowsOut(ii,1:length(peakInds))=powerRMnoiseAvRM(ii,peakInds);

    if ismember(ii,plotRangeInds) & ~isempty(plotTime)
        close all
        f0=figure('Position',[200 500 1200 500],'DefaultAxesFontSize',12,'visible',showPlot);
        plot(velOut(ii,:),powerRMnoiseRawPlot(ii,:),'-b')
        hold on
        plot([velOut(ii,1),velOut(ii,end)],[noiseThreshAll(ii),noiseThreshAll(ii)],'-c','LineWidth',2);
        plot([velOut(ii,1),velOut(ii,end)],[meanNoiseAll(ii),meanNoiseAll(ii)],'-k','LineWidth',2);
        plot(velOut(ii,:),powerRMnoiseAv(ii,:),'-g','LineWidth',1.5);
        plot(velOut(ii,:),powerRMnoiseAvRM(ii,:),'-r','LineWidth',1.5);
        plot(velOut(ii,:),powerSmoothAll(ii,:),'-k','LineWidth',1.5);
        scatter(velOut(ii,peakInds),powerSmoothAll(ii,peakInds),50,'filled','MarkerFaceColor','k');
        xlim([velOut(ii,1),velOut(ii,end)]);
        title(num2str(ii))
        hold off
        if ~plotAll
            set(gcf,'PaperPositionMode','auto')
            print(f0,[figdir,'spectra/spectra_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
        end
    end
end

%% Waterfall plot

if ~isempty(plotTime)
    f1=figure('Position',[200 500 1200 885],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');
    t=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);
    hold on
    surf(1:sampleNum,1:size(specDB,1),specDB,'edgecolor','none');
    view(2)
    xlim([1,sampleNum]);
    clim([-60,-10]);
    title('Raw')
    box on
    ylimits=s1.YLim;
    for kk=1:length(plotRangeInds)
        plot([1,sampleNum],[plotRangeInds(kk),plotRangeInds(kk)],'-m','LineWidth',1.5);
    end
    s1.YLim=ylimits;
    s1.SortMethod='childorder';

    s2=nexttile(2);
    hold on
    surf(1:sampleNum,1:size(powerRMnoiseAvRM,1),powerRMnoiseAvRM,'edgecolor','none');
    view(2)
    xlim([1,sampleNum]);
    s2.YLim=ylimits;
    clim([-60,-10]);
    title('Smoothed without noise')
    box on
    for kk=1:length(plotRangeInds)
        plot([1,sampleNum],[plotRangeInds(kk),plotRangeInds(kk)],'-k','LineWidth',1.5);
    end
    s3.SortMethod='childorder';

    s3=nexttile(3);
    hold on
    surf(1:sampleNum,1:size(powerRMnoiseRaw,1),powerRMnoiseRaw,'edgecolor','none');
    view(2)
    xlim([1,sampleNum]);
    s3.YLim=ylimits;
    clim([-60,-10]);
    title('Raw withouth noise')
    box on
    for kk=1:length(plotRangeInds)
        plot([1,sampleNum],[plotRangeInds(kk),plotRangeInds(kk)],'-k','LineWidth',1.5);
    end
    for aa=1:size(loopInds2,1)
        kk=loopInds2(aa);
        scatter(peaksIndsAll(kk,1),kk,30,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
        scatter(peaksIndsAll(kk,2),kk,30,'filled','MarkerFaceColor',[0.6,0.6,0.6],'MarkerEdgeColor','k');
    end
    s3.SortMethod='childorder';
    colorbar
    if ~plotAll
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,'waterfall/waterfall_',datestr(plotTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0');
    end
end
end