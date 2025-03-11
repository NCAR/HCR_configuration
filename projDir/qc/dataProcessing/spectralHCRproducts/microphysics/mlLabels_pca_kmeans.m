function [labelsInNum,kCenters,numLabel]=mlLabels_pca_kmeans(dataForML,numLabel,numComp,nRow,nCol)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% Principal component analysis PCA
disp('pca');
[~,~,~,~,explained]=pca(dataGood);
pareto(explained);
if isempty(numComp)
    % Find number of components
    numComp=input("Enter number of PCA components:");

    title(['Number of PCA components used: ',num2str(numComp)]);
end
print(['numComponents.png'],'-dpng','-r0');

[wcoeff,score,latent,tsquared,explained]=pca(dataGood,'NumComponents',numComp);

% K means
% Find Number of labels
if isempty(numLabel)
    minTestLabel=2;
    maxTestLabel=25;
    disp('Investigating number of k-means labels ...')
    sumdAll=nan(maxTestLabel-minTestLabel+1,1);
    for ii=minTestLabel:maxTestLabel
        [~,~,sumd]=kmeans(score,ii,'MaxIter',1000);
        sumdAll(ii-1)=sum(sumd);
    end

    f1 = figure('Position',[200 500 600 800],'DefaultAxesFontSize',12);
    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);
    scatter(1:length(sumdAll),sumdAll)
    s2=nexttile(2);
    plot(-diff(sumdAll))

    numLabel=input("Enter number of k-means labels:");

    title(['Number of k-means labels used: ',num2str(numLabel)]);

    print(f1,['numLabels.png'],'-dpng','-r0');
end

% Get k-means labels
disp('k-means');

[kInds,kCenters]=kmeans(score,numLabel,'MaxIter',1000);

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

end