function [err,resid]=noisePeaks_smoothingTest(specDB,velIn,data,widthC,aircVel,err,resid,figdir,plotTime)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries

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

if sampleNum~=987
    return
end

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

meanNoiseAll=nan(size(specDB,1),1);
noiseThreshAll=meanNoiseAll;

%% Remove noise
% Moving average
meanOverPoints=3; % Average over this number of points
secondMean=round(sampleNum/3);
movAv=movmedian(powerSpecLarge,meanOverPoints,2);
movAv2=movmedian(movAv,secondMean,2);

movAv2(:,1:round(sampleNum/3))=nan;
movAv2(:,end-round(sampleNum/3):end)=nan;

loopInds=find(any(~isnan(specDB),2));

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    thisMov=powerSpecLarge(ii,:);

    % Get one whole spectrum
    [~,minIndTest]=min(movAv2(ii,:));

    testPow=thisMov(minIndTest:minIndTest+sampleNum-1);
    testVel=velSpecLarge(minIndTest:minIndTest+sampleNum-1);

    % Find noise floor and noise threshold
    [noiseThreshAll(ii),meanNoiseAll(ii)]=findNoiseThresh(testPow,meanOverPoints);

    % Mean velocity
    noiseLinV=10.^(data.noise_v./10);
    sigInLin=10.^(testPow./10)-noiseLinV;

    % VEL
    meanVel=sum(sigInLin.*testVel,'omitmissing')/sum(sigInLin,'omitmissing');

    %filterAt=8;
    filterAt=round(0.00022396.*aircVel.^2-0.10542.*aircVel+18.132);
    filterAt=fillmissing(filterAt,'nearest');

    % Correct for aircraft width
    [err,errCat,sigWidthCorr,sigFiltered,signalIn1,signalIn2,sigFiltered1,sigFiltered2,inds1,inds2]= ...
        smoothingTest(filterAt,testPow,meanVel,widthC,testVel,sampleNum,err);

    % Calculate standar deviation of noise
    % for kk=1:length(filterAt)
    % residAdd=testPow-sigFiltered(filterAt(kk)-1,:);
    % end
    % resid=cat(2,resid,residAdd);
    resid=nan;


    if ismember(ii,plotRangeInds) & ~isempty(plotTime)

        errMean=mean(errCat,2);

        [~,errMinInd]=min(errMean);

        close all

        %close all
        f1=figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);
        t = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
        s1=nexttile(1);

        grid on
        box on
        hold on
        %plot(xVel,signalIn,'-b','LineWidth',2)
        plot(testVel(inds1),signalIn1,'-c','LineWidth',1)
        plot(testVel(inds2),signalIn2,'-k','LineWidth',1)
        xlim([testVel(1),testVel(end)]);

        legend('Split signal 1','Split signal 2')

        s2=nexttile(2);

        hold on
        plot(testVel,testPow,'-b','LineWidth',1)
        plot(testVel(inds1),sigFiltered1(filterAt-1,:),'-c','LineWidth',2)
        plot(testVel(inds2),sigFiltered2(filterAt-1,:),'-k','LineWidth',2)
        plot(testVel(inds1),sigFiltered1(errMinInd,:),'-y','LineWidth',1)
        plot(testVel(inds2),sigFiltered2(errMinInd,:),'-m','LineWidth',1)
        xlim([testVel(1),testVel(end)]);

        legend('Original signal','Filtered split 1','Filtered split 2', ...
            ['Filtered split 1 at nz ',num2str(errMinInd+1)],['Filtered split 2 at nz ',num2str(errMinInd+1)])

        grid on
        box on

        s3=nexttile(3);

        hold on
        plot(testVel,testPow,'-b','LineWidth',1)
        plot(testVel,sigFiltered(filterAt-1,:),'-g','LineWidth',2)
        plot(testVel,sigWidthCorr,'-r','LineWidth',2)
        hold off

        xlim([testVel(1),testVel(end)]);

        legend('Original signal','Filtered',['Filtered width corrected (',num2str(filterAt),' nz)']);

        grid on
        box on

        if ~plotAll
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,'spectra/zero',num2str(filterAt),'_everyOther_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
        end
    end
end
end