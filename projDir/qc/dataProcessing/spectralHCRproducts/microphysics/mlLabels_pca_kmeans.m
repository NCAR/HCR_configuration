function [labelsInNum,kCenters]=mlLabels_pca_kmeans(dataForML,numLabel,nRow,nCol)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% Principal component analysis
disp('pca');
% Find number of components
%[wcoeff2,score2,latent,tsquared,explained]=pca(dataGood);
%pareto(explained)

[wcoeff,score,latent,tsquared,explained]=pca(dataGood,'NumComponents',3);

% k-means labels
disp('k-means');
% Find Number of clusters
sumdAll=nan(numLabel-1,1);
for ii=2:numLabel
    [kInds,kCenters,sumd]=kmeans(score,ii,'MaxIter',1000);
    sumdAll(ii-1)=sum(sumd);
end
[kInds,kCenters]=kmeans(score,numLabel,'MaxIter',1000);

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

end