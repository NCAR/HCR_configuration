% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

numLabel=20;
figdirver='momsSpecAllT';

% Variables
%vars={'VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS','EDGE_EDGE_WIDTH','LEFT_SLOPE','RIGHT_SLOPE'}; %spec7vars
%vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS'}; %basicMoms
%vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS','TEMP'}; %basicMomsT
%vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS','LEFT_SLOPE','RIGHT_SLOPE','TEMP'}; %momsSpecT
%vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS','LEFT_SLOPE','RIGHT_SLOPE'}; %momsSpec
%vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS', ...
%    'LEFT_SLOPE','RIGHT_SLOPE','EDGE_EDGE_WIDTH','LEFT_EDGE_VEL','RIGHT_EDGE_VEL'}; %momsSpecAll
vars={'DBZ_MASKED','VEL_MASKED','WIDTH_SPEC','SKEWNESS','KURTOSIS', ...
    'LEFT_SLOPE','RIGHT_SLOPE','EDGE_EDGE_WIDTH','LEFT_EDGE_VEL','RIGHT_EDGE_VEL','TEMP'}; %momsSpecAllT

showPlot='off';

freqData='10hz_spec';

figdir=['/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full_spec/microphysics/train/', ...
    figdirver,'_',num2str(numLabel),'labels/'];

if ~exist(figdir,'dir')
    mkdir(figdir);
end

caseList=readtable('~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/trainMicroPhys.txt','Delimiter',' ');

startTimeAll=datetime(caseList{:,2:7});
endTimeAll=datetime(caseList{:,8:13});

%% Read data

disp(datetime)
disp('Aggregating data ...')

% Loop through files and aggregate data
caseStartInd=1;
for aa=1:size(caseList,1)
%for    aa=1:2

    startTime=startTimeAll(aa);
    endTime=endTimeAll(aa);

    if strcmp(caseList.Var1(aa),'sp')
        project='spicule'; %socrates, aristo, cset, otrec
        quality='qc1';
        qcVersion='v1.2';
    elseif strcmp(caseList.Var1(aa),'ot')
        project='otrec'; %socrates, aristo, cset, otrec
        quality='qc3';
        qcVersion='v3.2';
    elseif strcmp(caseList.Var1(aa),'so')
        project='socrates'; %socrates, aristo, cset, otrec
        quality='qc3';
        qcVersion='v3.2';
    elseif strcmp(caseList.Var1(aa),'cs')
        project='cset'; %socrates, aristo, cset, otrec
        quality='qc3';
        qcVersion='v3.1';
    end

    dataDir=HCRdir(project,quality,qcVersion,freqData);
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    for ii=1:length(vars)
        data.(vars{ii})=[];
    end
            
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
   
    if aa==1
        allVars=fields(data);
        specDataAll=data;
    else
        caseStartInd=cat(1,caseStartInd,length(specDataAll.time)+1);
        for ii=1:length(allVars)
            specDataAll.(allVars{ii})=cat(2,specDataAll.(allVars{ii}),data.(allVars{ii}));
        end
    end
end

%% Scale and prepare input

[dataForML,lims]=prepForML(specDataAll,vars);

%% K-means clustering

disp('Clustering ...');
[nRow,nCol]=size(specDataAll.(vars{1}));
%[labelsKmeans,centersKmeans,centKmReScaled]=mlLabels(dataForML,numLabel,nRow,nCol,lims,vars);
[labelsKmeans,centersKmeans,centKmReScaled]=mlLabels_sil(dataForML,numLabel,nRow,nCol,lims,vars);

save([figdir,'centers.mat'],'centersKmeans','vars');

%% Set up analysis

uTypes=unique(caseList.Var14(:));

for ii=1:length(uTypes)
    hists.(uTypes{ii})=[];
end

histAll=[];

%% Plot
caseEndInd=caseStartInd-1;
caseEndInd(1)=[];
caseEndInd=cat(1,caseEndInd,length(specDataAll.time));

edges=0.5:1:numLabel+0.5;
%cmap=cat(1,[0,0,0],tab20(numLabel));
cmap=cat(1,[0,0,0],distinguishable_colors(numLabel,'k'));

disp('Plotting ...');
for aa=1:size(caseList,1)

    close all

    for jj=1:length(vars);
        specDataPlot.(vars{jj})=specDataAll.(vars{jj})(:,caseStartInd(aa):caseEndInd(aa));
    end
    specDataPlot.time=specDataAll.time(caseStartInd(aa):caseEndInd(aa));
    specDataPlot.asl=specDataAll.asl(:,caseStartInd(aa):caseEndInd(aa));
    specDataPlot.labelsKmeans=labelsKmeans(:,caseStartInd(aa):caseEndInd(aa));

    % Histogram
    histC=histcounts(specDataPlot.labelsKmeans(~isnan(specDataPlot.labelsKmeans)),edges);
    histAll=cat(1,histAll,histC);
    histC=histC./sum(histC).*100;
    hists.(caseList.Var14{aa})=cat(1,hists.(caseList.Var14{aa}),histC);
            
    aslGood=specDataPlot.asl(~isnan(specDataPlot.labelsKmeans))./1000;
    ylims=[0,max(aslGood)+0.5];

    specDataPlot.labelsKmeans(isnan(specDataPlot.labelsKmeans))=-99;

    f1 = figure('Position',[200 500 1200 1150],'DefaultAxesFontSize',12,'visible',showPlot);

    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    hold on
    surf(specDataPlot.time,specDataPlot.asl./1000,specDataPlot.labelsKmeans,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([-0.5,numLabel+0.5]);
    s1.Colormap=cmap;
    colorbar
    grid on
    box on
    title('k-means labels')
    ylim(ylims);
    xlim([specDataPlot.time(1),specDataPlot.time(end)]);

    s2=nexttile(2);
    b=bar(1:numLabel,histC,'FaceColor','flat');
    for kk = 1:numLabel
        b.CData(kk,:)=cmap(kk+1,:);
    end

    xlabel('Label');
    ylabel('Percent (%)');

    xlim([0.5,numLabel+0.5]);
    xticks(1:numLabel);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'labels_',caseList.Var14{aa},'_',caseList.Var1{aa},'_', ...
        datestr(specDataPlot.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(specDataPlot.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
    
    % Plot input fields
    plotFields(specDataPlot,vars,caseList.Var14{aa},caseList.Var1{aa},figdir,ylims,showPlot);
    
end

%% Histograms

f1 = figure('Position',[200 500 1400 1200],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(4,1,'TileSpacing','tight','Padding','compact');

colT=[1,0,0;0,0,0.5;0,0,1;0,0.5,0;0,1,0;0.2,0.2,0.2];
addX=-0.25:0.1:0.25;
x=1:numLabel;
x=x+addX';

s1=nexttile(1);
hold on

for ii=1:length(uTypes)
    thisType=hists.(uTypes{ii});
    thisType(thisType>40)=40;
    for jj=1:size(thisType,1)
        scatter(x(ii,:),thisType(jj,:),'filled','MarkerFaceColor',colT(ii,:));
    end
end

xlim([0,numLabel+1]);
xticks(1:numLabel);
xlabel('Label');
ylabel('Percent (%)');
box on


s2=nexttile(2);
histCA=sum(histAll,1);
b=bar(1:numLabel,histCA./sum(histCA).*100,'FaceColor','flat');
for kk = 1:numLabel
    b.CData(kk,:)=cmap(kk+1,:);
end

xlabel('Label');
xticks(1:numLabel);
ylabel('Percent (%)');

xlim([0,numLabel+1]);

box on

s3=nexttile(3,[2 1]);
hold on
scatter(1:numLabel,centersKmeans',60,'*');
plot([0,numLabel+1],[1,1],'-k')

xlim([0,numLabel+1]);
ylim([0,1.08])
box on
legend(vars,'Interpreter','none','Location','north','Orientation','horizontal')
%colororder("gem12")
colororder(cmap(1:length(vars),:));

xticks(1:numLabel);
xlabel('Label');
yticks(0:0.1:1);
ylabel('Normalized variable value')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'histogram.png'],'-dpng','-r0');
    