function powerRMnoiseOut=rmNoiseSpec_HS(powerAdj,vNoise,plotTime)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
powerRMnoiseOut=nan(size(powerAdj));

sampleNumR=size(powerAdj,2);

% Moving average
meanOverPoints=5; % Average over this number of points
movAv=movmedian(powerAdj,meanOverPoints,2);

startNoise=10^((vNoise+15)/10);

powerRMnoiseAv=nan(size(powerAdj));

plotRangeInds=[20:20:700];
loopInds=find(any(~isnan(powerAdj),2));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/figsMultiVel/cases/spectra/';
showPlot='off';
for aa=1:size(loopInds,1)
    ii=loopInds(aa);

    thisMov=movAv(ii,:);

    % Find noise floor and noise threshold
    thisMovLin=10.^(thisMov./10);
    [noiseThreshM,meanNoiseM]=findNoiseThresh(thisMovLin,meanOverPoints,startNoise);

    % Remove noise
    thisMovRM=thisMov;
    thisMovRM(thisMovRM<noiseThreshM)=nan;
    thisMask1=~isnan(thisMovRM);
    lastwarn(''); % Clear last warning message
    thisMask=bwareafilt(thisMask1,1);

    % bwareafilt creates a warning if more than one region is the largest.
    % Check if a warning was created, if yes, find the 5 largest and pick
    % the one with the highest power
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        warnStruct=warning('off',warnId);
        
        thisMask=bwareafilt(thisMask1,5);
        thisMovRM(thisMask==0)=nan;

        [~,maxInd]=max(thisMovRM,[],2,'omitmissing');
        specRegs=bwconncomp(thisMask);

        for jj=1:specRegs.NumObjects
            if ~ismember(maxInd,specRegs.PixelIdxList{jj})
                thisMask(specRegs.PixelIdxList{jj})=0;
            end
        end
    end

    thisMovRM(thisMask==0)=nan;
    powerRMnoiseAv(ii,:)=thisMovRM;
    powerRMnoiseOut(ii,thisMask==1)=powerAdj(ii,thisMask==1);

    if ismember(ii,plotRangeInds) & ~isempty(plotTime)
        close all
        f0=figure('visible',showPlot);
        plot(1:sampleNumR,powerAdj(ii,:),'-b')
        hold on
        plot(1:sampleNumR,thisMov,'-g','LineWidth',1.5);
        plot(1:sampleNumR,thisMovRM,'-r','LineWidth',1.5);
        plot([1,sampleNumR],[noiseThreshM,noiseThreshM],'-c','LineWidth',2);
        plot([1,sampleNumR],[meanNoiseM,meanNoiseM],'-k','LineWidth',2);
        xlim([1,sampleNumR]);
        title(num2str(ii))
        hold off
        set(gcf,'PaperPositionMode','auto')
        print(f0,[figdir,'spectra_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
    end
end

if ~isempty(plotTime)
    f1=figure('Position',[200 500 1200 885],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');
    t=tiledlayout(1,3,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);
    hold on
    surf(1:sampleNumR,1:size(powerAdj,1),powerAdj,'edgecolor','none');
    view(2)
    xlim([1,sampleNumR]);
    clim([-60,-10]);
    title('Raw')
    box on
    ylimits=s1.YLim;
    for kk=1:length(plotRangeInds)
        plot([1,sampleNumR],[plotRangeInds(kk),plotRangeInds(kk)],'-m','LineWidth',1.5);
    end
    s1.YLim=ylimits;
    s1.SortMethod='childorder';

    s2=nexttile(2);
    hold on
    surf(1:sampleNumR,1:size(powerAdj,1),powerRMnoiseAv,'edgecolor','none');
    view(2)
    xlim([1,sampleNumR]);
    s2.YLim=ylimits;
    clim([-60,-10]);
    title('Smoothed without noise')
    box on
    for kk=1:length(plotRangeInds)
        plot([1,sampleNumR],[plotRangeInds(kk),plotRangeInds(kk)],'-k','LineWidth',1.5);
    end
    s3.SortMethod='childorder';

    s3=nexttile(3);
    hold on
    surf(1:sampleNumR,1:size(powerAdj,1),powerRMnoiseOut,'edgecolor','none');
    view(2)
    xlim([1,sampleNumR]);
    s3.YLim=ylimits;
    clim([-60,-10]);
    title('Raw withouth noise')
    box on
    for kk=1:length(plotRangeInds)
        plot([1,sampleNumR],[plotRangeInds(kk),plotRangeInds(kk)],'-k','LineWidth',1.5);
    end
    s3.SortMethod='childorder';
    colorbar

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'waterfall_',datestr(plotTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0');
end
end