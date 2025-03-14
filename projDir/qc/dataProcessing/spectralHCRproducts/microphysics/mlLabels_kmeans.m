function [labelsInNum,kCenters,kCentReScaled,numLabel]=mlLabels_kmeans(dataForML,numLabel,nRow,nCol,lims,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

%% Correlation
C=corr(dataGood,dataGood);
f1 = figure('Position',[200 500 900 800],'DefaultAxesFontSize',12);
h=heatmap(C,'XData',vars,'YData',vars,'Interpreter','none','CellLabelFormat','%.2f');
colormap(turbo(20));
clim([-1,1]);
title('Correlation')
print('correlation.png','-dpng','-r0');

%% Principal component analysis PCA
disp('pca');
[coeff,~,~,~,explained]=pca(dataGood);

f1 = figure('Position',[200 500 900 1000],'DefaultAxesFontSize',12);
t = tiledlayout(3,1,'TileSpacing','compact','Padding','tight');
s1=nexttile(1);
pareto(explained,1);
xlabel('Principal components')
s2=nexttile(2,[2,1]);
heatmap(coeff,'YData',vars,'Interpreter','none','CellLabelFormat','%.2f')
colormap(turbo(length(vars)));
clim([-1,1]);
xlabel('Principal components')
title('Principal component coefficients')
print('pcaCoeffs.png','-dpng','-r0');

%% Find Number of k-means labels
if isempty(numLabel)
    minTestLabel=2;
    maxTestLabel=25;
    disp('Investigating number of k-means labels ...')

    % Subsample for performance
    dataGoodSub=dataGood(1:100:size(dataGood,1),:);
    % Elbow method
    disp('Ellbow');
    evaEll=nan(maxTestLabel-minTestLabel+1,1);
    for ii=minTestLabel:maxTestLabel
        [~,~,sumd]=kmeans(dataGoodSub,ii,'Replicates',5,'MaxIter',1000);
        evaEll(ii-1)=sum(sumd);
    end

    % With evalclusters function
    clustFcn=@(X,K) kmeans(X, K,'Replicates',5,'MaxIter',1000);
    disp('Silhouette');
    evaSil=evalclusters(dataGoodSub,clustFcn,'silhouette','KList',minTestLabel:maxTestLabel);
    disp('Davies Bouldin');
    evaDB=evalclusters(dataGoodSub,clustFcn,'DaviesBouldin','KList',minTestLabel:maxTestLabel);
    disp('Gap');
    dataGoodSubG=dataGood(1:5000:size(dataGood,1),:);
    evaGap=evalclusters(dataGoodSubG,clustFcn,'gap','KList',minTestLabel:maxTestLabel);

    f1 = figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12);
    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);
    hold on
    plot(minTestLabel:maxTestLabel,rescale(evaEll,-1,1),'b-o');
    plot(minTestLabel:maxTestLabel,rescale(evaSil.CriterionValues,-1,1),'r-x');
    plot(minTestLabel:maxTestLabel,rescale(evaDB.CriterionValues,-1,1),'g-*');
    plot(minTestLabel:maxTestLabel,rescale(evaGap.CriterionValues,-1,1),'k-s');

    xticks(1:maxTestLabel);
    grid on
    box on

    legend('Ellbow (diff)','Silhouette (max)','Davies Bouldin (min)','Gap (max)', ...
        'Location','northoutside','Orientation','horizontal')
    
    s2=nexttile(2);
    hold on
    plot(minTestLabel+0.5:maxTestLabel-0.5,rescale(diff(evaEll),-1,1),'b-o');
    plot(minTestLabel+0.5:maxTestLabel-0.5,rescale(-diff(evaSil.CriterionValues),-1,1),'r-x');
    plot(minTestLabel+0.5:maxTestLabel-0.5,rescale(diff(evaDB.CriterionValues),-1,1),'g-*');
    plot(minTestLabel+0.5:maxTestLabel-0.5,rescale(diff(evaGap.CriterionValues),-1,1),'k-s');

    xticks(1:maxTestLabel);
    grid on
    box on

    title('Difference')

    % Get labels from user
    numLabel=input("Enter number of k-means labels:");

    title(s1,['Number of k-means labels used: ',num2str(numLabel)]);

    print(f1,['numLabels.png'],'-dpng','-r0');
end

% k-means labels
disp('k-means');
[kInds,kCenters]=kmeans(dataGood,numLabel,'MaxIter',1000);

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

% Re-scale centers
kCentReScaled=nan(size(kCenters));
for ii=1:length(vars)
    kCentReScaled(:,ii)=kCenters(:,ii).*(lims.(vars{ii})(2)-lims.(vars{ii})(1))+lims.(vars{ii})(1);
    %kCentReScaled2(:,ii)=rescale(kCenters(:,ii),lims.(vars{ii})(1),lims.(vars{ii})(2));
end
end