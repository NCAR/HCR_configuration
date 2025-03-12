function [labelsInNum,kCenters,numLabel]=mlLabels_pca_kmeans(dataForML,numLabel,numComp,nRow,nCol,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% Correlation
C=corr(dataGood,dataGood);
f1 = figure('Position',[200 500 900 800],'DefaultAxesFontSize',12);
heatmap(C,'XData',vars,'YData',vars,'Interpreter','none','CellLabelFormat','%.2f')
colormap(turbo(20));
clim([-1,1]);
title('Correlation')
print(['correlation.png'],'-dpng','-r0');

% Principal component analysis PCA
disp('pca');
[coeff,~,~,~,explained]=pca(dataGood);
pareto(explained,1);
if isempty(numComp)
    % Find number of components
    numComp=input("Enter number of PCA components:");
end
title(['Number of PCA components used: ',num2str(numComp)]);
print(['numComponents.png'],'-dpng','-r0');

f1 = figure('Position',[200 500 900 800],'DefaultAxesFontSize',12);
heatmap(coeff,'YData',vars,'Interpreter','none','CellLabelFormat','%.2f')
colormap(turbo(length(vars)));
clim([-1,1]);
xlabel('Principal components')
title('Principal component coefficients')
print(['pcaCoeffs.png'],'-dpng','-r0');

[coeff,score,latent,tsquared,explained]=pca(dataGood,'NumComponents',numComp);

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
    scatter(minTestLabel:maxTestLabel,sumdAll);
    s2=nexttile(2);
    scatter(minTestLabel+0.5:maxTestLabel-0.5,-diff(sumdAll));

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