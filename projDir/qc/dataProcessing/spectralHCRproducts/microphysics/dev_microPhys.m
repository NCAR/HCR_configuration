% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.2';

dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

indir=[dataDirCF(1:end-5),'microPhys/matFiles/'];
figdir=[dataDirCF(1:end-5),'microPhys/training/'];

showPlot='on';

% casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/airMotion_',project,'.txt'];
fileList=dir([indir,'*.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through filed

for aa=1:length(fileList)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(fileList))]);
   
    load([fileList(aa).folder,'/',fileList(aa).name]);

    %% Scale and prepare input
    vars={'DBZ','VEL','WIDTH','SKEW','KURT','MELTING_LAYER'};
    dataForML=prepForML(specData,vars);

    %% K-means clustering

    numLabel=10;
    [nRow,nCol]=size(specData.(vars{1}));
    [labelsInNum,labelsOptEll,labNumEll]=mlLabels(dataForML,numLabel,nRow,nCol);

    %% Plot labels

    close all

    disp('Plotting ...');

    aslGood=specData.asl(~isnan(specData.DBZ))./1000;
    ylims=[0,max(aslGood)+0.5];

    f1 = figure('Position',[200 500 700 800],'DefaultAxesFontSize',12,'visible',showPlot);

    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    hold on
    surf(specData.time,specData.asl./1000,labelsInNum,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0.5,numLabel+0.5]);
    s1.Colormap=turbo(numLabel);
    colorbar
    grid on
    box on
    title('k-means labels')
    ylim(ylims);
    xlim([specData.time(1),specData.time(end)]);

    s2=nexttile(2);

    hold on
    surf(specData.time,specData.asl./1000,labelsOptEll,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0.5,labNumEll+0.5]);
    s2.Colormap=turbo(labNumEll);
    colorbar
    grid on
    box on
    title('k-means labels ellbow optimized')
    ylim(ylims);
    xlim([specData.time(1),specData.time(end)]);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_kmLabels_',datestr(specData.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(specData.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    
    %% Plot moments and spec params
    %plotAllMoments(dataCF,momentsSpecSmoothCorr,figdir,project,showPlot);
    %plotSpecParams(momentsSpecSmoothCorr,figdir,project,showPlot);

end
