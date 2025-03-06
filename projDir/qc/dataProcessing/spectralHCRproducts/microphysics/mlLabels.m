function [labelsInNum,kCenters,kCentReScaled]=mlLabels(dataForML,numLabel,nRow,nCol,lims,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% k-means labels
disp('k-means');
[kInds,kCenters]=kmeans(dataGood,numLabel,'MaxIter',1000);

% distancesFromOrigin = sqrt(kCenters(:, 1) .^ 2 + kCenters(:, 2) .^2);
% 
% [~,sortOrder] = sort(distancesFromOrigin, 'ascend'); % Sort x values of centroids.
% kIndsS = zeros(size(dataGood,1), 1);
% % For each class, find out where it is
% for k = 1 : size(kCenters, 1)
%     currentClassLocations = kInds == k;
%     newClassNumber = find(k == sortOrder);	% Find index in sortOrder where this class number appears.
%     kIndsS(currentClassLocations) = newClassNumber;
% end

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

% Re-scale centers
kCentReScaled=nan(size(kCenters));
for ii=1:length(vars)
    kCentReScaled(:,ii)=kCenters(:,ii).*(lims.(vars{ii})(2)-lims.(vars{ii})(1))+lims.(vars{ii})(1);
end
end