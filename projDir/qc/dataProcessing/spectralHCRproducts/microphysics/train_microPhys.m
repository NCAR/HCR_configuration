% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='ts'; %field, qc1, or qc2
qualityCF='qc1';
freqData='10hz';
qcVersion='v1.2';

numFiles=30; % Maximum 25

dataDirCF=HCRdir(project,qualityCF,qcVersion,freqData);

%indir=[dataDirCF(1:end-5),'microPhys/matFiles/'];
figdir=[dataDirCF(1:end-5),'microPhys/training/'];

showPlot='on';

fileList=readtable('~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/trainMicroPhys.txt','Delimiter',' ');
%fileList=dir([indir,'*.mat']);

gapL=round(size(fileList,1)/numFiles);
fileInds=1:gapL:size(fileList,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Aggregating data ...')

% Loop through files and aggregate data
caseStartInd=1;
for bb=1:length(fileInds)

    aa=fileInds(bb);

    %load([fileList(aa).folder,'/',fileList(aa).name]);
    load(fileList.File{aa});
    
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
%vars={'DBZ','VEL','WIDTH','SKEW','KURT'}; % Version 1 (10), Version 4 (15), Version 5 (20)
%vars={'DBZ','VEL','WIDTH','SKEW','KURT','EE_WIDTH','L_SLOPE','R_SLOPE','MELTING_LAYER'}; % Version 2 (20)
%vars={'DBZ','VEL','WIDTH','SKEW','KURT','EE_WIDTH','L_SLOPE','R_SLOPE'}; % Version 3 (15)
%vars={'VEL','WIDTH','SKEW','KURT'}; % Version 6 (20),
vars={'VEL','WIDTH','SKEW','KURT','EE_WIDTH','L_SLOPE','R_SLOPE'}; % Version 7 (20)

dataForML=prepForML(specDataAll,vars);

%% K-means clustering

disp('Clustering ...');
numLabel=20;
[nRow,nCol]=size(specDataAll.(vars{1}));
[labelsKmeans,centersKmeans]=mlLabels(dataForML,numLabel,nRow,nCol);

save([figdir,'centers.mat'],'centersKmeans','vars');

%% Plot
caseEndInd=caseStartInd-1;
caseEndInd(1)=[];
caseEndInd=cat(1,caseEndInd,length(specDataAll.time));

disp('Plotting ...');
for aa=1:length(fileInds)

    close all

    for jj=1:length(vars);
        specDataPlot.(vars{jj})=specDataAll.(vars{jj})(:,caseStartInd(aa):caseEndInd(aa));
    end
    specDataPlot.time=specDataAll.time(caseStartInd(aa):caseEndInd(aa));
    specDataPlot.asl=specDataAll.asl(:,caseStartInd(aa):caseEndInd(aa));
    specDataPlot.labelsKmeans=labelsKmeans(:,caseStartInd(aa):caseEndInd(aa));
        
    aslGood=specDataPlot.asl(~isnan(specDataPlot.labelsKmeans))./1000;
    ylims=[0,max(aslGood)+0.5];

    f1 = figure('Position',[200 500 700 300],'DefaultAxesFontSize',12,'visible',showPlot);

    t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    hold on
    surf(specDataPlot.time,specDataPlot.asl./1000,specDataPlot.labelsKmeans,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([0.5,numLabel+0.5]);
    s1.Colormap=tab20(numLabel);
    colorbar
    grid on
    box on
    title('k-means labels')
    ylim(ylims);
    xlim([specDataPlot.time(1),specDataPlot.time(end)]);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'labels_',datestr(specDataPlot.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(specDataPlot.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    
    % Plot input fields
    plotFields(specDataPlot,vars,figdir,ylims,showPlot);
    
end
