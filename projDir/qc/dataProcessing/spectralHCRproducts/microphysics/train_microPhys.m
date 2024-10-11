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

disp('Aggregating data ...')

% Loop through files and aggregate data
caseStartInd=1;
for aa=1:length(fileList)

    load([fileList(aa).folder,'/',fileList(aa).name]);

    if aa==1
        allVars=fields(specData);
        specDataAll=specData;
    else
        caseStartInd=cat(1,caseStartInd,length(specDataAll.time)+1);
        for ii=1:length(allVars)
            specDataAll.(allVars{ii})=cat(2,specDataAll.(allVars{ii}),specData.(allVars{ii}));
        end
    end
end

%% Scale and prepare input
vars={'DBZ','VEL','WIDTH','SKEW','KURT'};
dataForML=prepForML(specDataAll,vars);

%% K-means clustering

disp('Clustering ...');
numLabel=10;
[nRow,nCol]=size(specDataAll.(vars{1}));
[labelsInNum,labelsOptEll,labNumEll]=mlLabels(dataForML,numLabel,nRow,nCol);

%% Plot labels
caseEndInd=caseStartInd-1;
caseEndInd(1)=[];
caseEndInd=cat(1,caseEndInd,length(specDataAll.time));

disp('Plotting ...');
for aa=1:length(fileList)
    close all

    time=specDataAll.time(caseStartInd(aa):caseEndInd(aa));
    asl=specDataAll.asl(:,caseStartInd(aa):caseEndInd(aa));
    labelsInNumC=labelsInNum(:,caseStartInd(aa):caseEndInd(aa));
    labelsOptEllC=labelsOptEll(:,caseStartInd(aa):caseEndInd(aa));
    
    aslGood=asl(~isnan(labelsInNumC))./1000;
    ylims=[0,max(aslGood)+0.5];

    f1 = figure('Position',[200 500 700 800],'DefaultAxesFontSize',12,'visible',showPlot);

    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    hold on
    surf(time,asl./1000,labelsInNumC,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0.5,numLabel+0.5]);
    s1.Colormap=turbo(numLabel);
    colorbar
    grid on
    box on
    title('k-means labels')
    ylim(ylims);
    xlim([time(1),time(end)]);

    s2=nexttile(2);

    hold on
    surf(time,asl./1000,labelsOptEllC,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0.5,labNumEll+0.5]);
    s2.Colormap=turbo(labNumEll);
    colorbar
    grid on
    box on
    title('k-means labels ellbow optimized')
    ylim(ylims);
    xlim([time(1),time(end)]);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_kmLabels_',datestr(time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    
    %% Plot moments and spec params
    %plotAllMoments(dataCF,momentsSpecSmoothCorr,figdir,project,showPlot);
    %plotSpecParams(momentsSpecSmoothCorr,figdir,project,showPlot);

end
