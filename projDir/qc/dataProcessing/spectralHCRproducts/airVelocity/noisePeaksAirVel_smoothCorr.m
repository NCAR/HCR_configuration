function [powerRMnoiseAvRM,powerRMnoise,powerRMnoiseAvRMS,velOut,velOutS]=noisePeaksAirVel_smoothCorr(specDB,velIn,data,widthC,figdir,plotTime)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
% Adjust spectra so they fit in the boundaries
powerRMnoiseAvRM=[];
powerRMnoise=[];
powerRMnoiseAvRMS=[];
velOut=[];
velOutS=[];


% Decide if and what to plot
plotAll=1; % Set to 1 if everything should be plotted. Plots won't be saved.
showPlot='on';

if plotAll
    plotRangeInds=18:10:size(specDB,1);
    plotTime=1;
else
    plotRangeInds=20:20:size(specDB,1);
end

sampleNum=length(data.time);

duplicateSpec=7;

% Add spectra side by side
powerSpecLarge=repmat(specDB,1,duplicateSpec);

velSpecLarge=-duplicateSpec*pi:2*pi/(sampleNum):duplicateSpec*pi;
velSpecLarge=velSpecLarge(1:end-1).*data.lambda./(4*pi.*repmat(data.prt,1,duplicateSpec));

meanNoiseAll=nan(size(specDB,1),1);
noiseThreshAll=meanNoiseAll;

%% Remove noise
% Moving average
meanOverPoints=round(sampleNum/3);
movAv=movmedian(powerSpecLarge,meanOverPoints,2);

movAv(:,1:round(sampleNum/3))=nan;
movAv(:,end-round(sampleNum/3):end)=nan;

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

noiseStd=5.52;

for aa=1:size(loopInds,1)
    ii=loopInds(aa); % ii is the range index

    % Find noise floor and remove noise
    [sigWidthCorrRMnoise,noiseThresh]=noiseFloorRM(sigWidthCorr(ii,:),noiseStd,sigPeaks(ii,:),sigValleys(ii,:),testVel(ii,:),fakeMeanVel(ii));

    if ismember(ii,plotRangeInds) & ~isempty(plotTime)

        close all

        %close all
        f1=figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);
        t = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');
       
        s1=nexttile(1);

        hold on
        plot(testVel(ii,:),testPow(ii,:),'-b','LineWidth',1)
        plot(testVel(ii,:),sigFiltered(ii,:),'-g','LineWidth',2)
        plot(testVel(ii,:),sigWidthCorr(ii,:),'-k','LineWidth',2)
        plot(testVel(ii,:),sigWidthCorrRMnoise,'-r','LineWidth',2)
        hold off

        xlim([testVel(ii,1),testVel(ii,end)]);

        legend({'Original signal','Filtered','Filtered width corrected','Without noise'},'Location','northoutside','Orientation','horizontal');

        grid on
        box on

        if ~plotAll
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,'spectra/spectra_',datestr(plotTime,'yyyymmdd_HHMMSS_'),num2str(ii),'.png'],'-dpng','-r0');
        end
    end
end
end