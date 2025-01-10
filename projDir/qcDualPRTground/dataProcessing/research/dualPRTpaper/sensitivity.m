clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

% Load other
load(['/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/sensitivity/SNRs.mat']);

%% Spatial coverage

countLong=sum(~isnan(longSNR(:)));
countShort=sum(~isnan(shortSNR(:)));
areaFracDiff=countLong/countShort;

disp(['Long pulse has ',num2str((areaFracDiff-1)*100),'% more data.']);

%% Difference
% Total
diffMat=longSNR-shortSNR;
totalValDiff=median(diffMat(:),'omitmissing');
totalStdDiff=std(diffMat(:),'omitmissing');

disp(['Median SNR difference: ',num2str(totalValDiff),' +/- ',num2str(totalStdDiff),' dB.']);

% By range
rangeValDiff=mean(diffMat,2,'omitmissing');
nonNan=sum(~isnan(diffMat),2);

%% Create data subsets
allLongSNR=longSNR(~isnan(longSNR));
allLongDBZ=longDBZ(~isnan(longDBZ));
allShortSNR=shortSNR(~isnan(shortSNR));
allShortDBZ=shortDBZ(~isnan(shortDBZ));

allLongSNRcens=longSNR(~isnan(shortSNR));
allShortSNRcens=shortSNR(~isnan(longSNR));

alts={'1','5','10'};

for ii=1:length(alts)
    altNum=str2num(alts{ii});
    [~,rangeAltInd]=min(abs(altNum*1000-range1));

    smallMat=longSNR(rangeAltInd-1:rangeAltInd+1,:);
    longSNRalt.(['a',alts{ii}])=smallMat(~isnan(smallMat));
    smallMat=shortSNR(rangeAltInd-1:rangeAltInd+1,:);
    shortSNRalt.(['a',alts{ii}])=smallMat(~isnan(smallMat));

    smallMat=longDBZ(rangeAltInd-1:rangeAltInd+1,:);
    longDBZalt.(['a',alts{ii}])=smallMat(~isnan(smallMat));
    smallMat=shortDBZ(rangeAltInd-1:rangeAltInd+1,:);
    shortDBZalt.(['a',alts{ii}])=smallMat(~isnan(smallMat));
end

%% Bin data

edgesIn=[-100:1:100];
edgesPlot=edgesIn(1:end-1)+(edgesIn(2)-edgesIn(1))/2;

allLongSNRbins=histcounts(allLongSNR,edgesIn);
allLongDBZbins=histcounts(allLongDBZ,edgesIn);
allShortSNRbins=histcounts(allShortSNR,edgesIn);
allShortDBZbins=histcounts(allShortDBZ,edgesIn);

allLongSNRbinsCens=histcounts(allLongSNRcens,edgesIn);
allShortSNRbinsCens=histcounts(allShortSNRcens,edgesIn);

for ii=1:length(alts)
    longSNRaltBin.(['a',alts{ii}])=histcounts(longSNRalt.(['a',alts{ii}]),edgesIn);
    longDBZaltBin.(['a',alts{ii}])=histcounts(longDBZalt.(['a',alts{ii}]),edgesIn);
    shortSNRaltBin.(['a',alts{ii}])=histcounts(shortSNRalt.(['a',alts{ii}]),edgesIn);
    shortDBZaltBin.(['a',alts{ii}])=histcounts(shortDBZalt.(['a',alts{ii}]),edgesIn);
end

%% Plot reflectivity
close all
f1 = figure('Position',[200 500 900 1000],'DefaultAxesFontSize',12,'renderer','painters');
t = tiledlayout(4,2,'TileSpacing','tight','Padding','compact');

s1=nexttile(1);
hold on
plot(edgesPlot,allLongDBZbins./sum(allLongDBZbins).*100,'-b','LineWidth',2);
plot(edgesPlot,allShortDBZbins/sum(allLongDBZbins).*100,'-g','LineWidth',2);

xlim([-55,20]);

legend('Long','Short','Location','northwest');

grid on
box on

xlabel('Reflectivity (dBZ)');
ylabel('Percent (%)')

title('(a) All reflectivities')

s2=nexttile(2);
hold on
plot(edgesPlot,allLongDBZbins./sum(allLongDBZbins).*100,'-b','LineWidth',2);
plot(edgesPlot,allShortDBZbins./sum(allLongDBZbins).*100,'-g','LineWidth',2);

xlim([-55,-40]);

legend('Long','Short','Location','northwest');

grid on
box on

xlabel('Reflectivity (dBZ)');
ylabel('Percent (%)');

title('(b) All reflectivities (zoomed)')

% Other panels
subplots={'(c)','(d)','(e)','(f)','(g)','(h)'};
xlims=[-50,20;
    -50,-35;
    -35,20;
    -35,-20;
    -30,15;
    -30,-15];

ylims=[0,5;
    0,3.9;
    0,6;
    0,2.9;
    0,8;
    0,2.2];

tileNum=3;

for ii=1:length(alts)
    
    s=nexttile(tileNum);
    hold on
    plot(edgesPlot,longDBZaltBin.(['a',alts{ii}])./sum(longDBZaltBin.(['a',alts{ii}])).*100,'-b','LineWidth',2);
    plot(edgesPlot,shortDBZaltBin.(['a',alts{ii}])./sum(longDBZaltBin.(['a',alts{ii}])).*100,'-g','LineWidth',2);

    xlim(xlims(tileNum-2,:));
    ylim(ylims(tileNum-2,:));

    %legend('Long','Short','Location','northwest');

    grid on
    box on

    xlabel('Reflectivity (dBZ)');
    ylabel('Percent (%)')

    title([subplots{tileNum-2},' Reflectivities at ',alts{ii},' km range']);

    tileNum=tileNum+1;

    s=nexttile(tileNum);
    hold on
    plot(edgesPlot,longDBZaltBin.(['a',alts{ii}])./sum(longDBZaltBin.(['a',alts{ii}])).*100,'-b','LineWidth',2);
    plot(edgesPlot,shortDBZaltBin.(['a',alts{ii}])./sum(longDBZaltBin.(['a',alts{ii}])).*100,'-g','LineWidth',2);

    xlim(xlims(tileNum-2,:));
    ylim(ylims(tileNum-2,:));

    %legend('Long','Short','Location','northwest');

    grid on
    box on

    xlabel('Reflectivity (dBZ)');
    ylabel('Percent (%)')

    title([subplots{tileNum-2},' Reflectivities at ',alts{ii},' km range (zoomed)']);

    tileNum=tileNum+1;
end

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'sensitivityDBZ.png'],'-dpng','-r0')

%% Plot SNR
close all
f1 = figure('Position',[200 500 600 400],'DefaultAxesFontSize',12,'renderer','painters');
t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact');

s1=nexttile(1);
hold on
plot(edgesPlot,allLongSNRbinsCens./sum(allLongSNRbinsCens).*100,'-b','LineWidth',2);
plot(edgesPlot,allShortSNRbinsCens/sum(allShortSNRbinsCens).*100,'-g','LineWidth',2);

xlim([-20,70]);

legend('Long','Short','Location','northeast');

grid on
box on

xlabel('SNR (dB)');
ylabel('Percent (%)')

title('(a) All SNRs')


set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'sensitivitySNR.png'],'-dpng','-r0')