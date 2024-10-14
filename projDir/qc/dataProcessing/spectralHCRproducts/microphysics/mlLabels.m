function [labelsInNum]=mlLabels(dataForML,numLabel,nRow,nCol)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% k-means labels
disp('k-means');
[kInds,kCenters]=kmeans(dataGood,numLabel,'MaxIter',1000);

distancesFromOrigin = sqrt(kCenters(:, 1) .^ 2 + kCenters(:, 2) .^2);

[~,sortOrder] = sort(distancesFromOrigin, 'ascend'); % Sort x values of centroids.
kIndsS = zeros(size(dataGood,1), 1);
% For each class, find out where it is
for k = 1 : size(kCenters, 1)
    currentClassLocations = kInds == k;
    newClassNumber = find(k == sortOrder);	% Find index in sortOrder where this class number appears.
    kIndsS(currentClassLocations) = newClassNumber;
end

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kIndsS;

labelsInNum=reshape(kIndsAll,nRow,nCol);
end