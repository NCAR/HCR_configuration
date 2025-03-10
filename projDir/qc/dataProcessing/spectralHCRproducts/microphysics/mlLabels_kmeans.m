function [labelsInNum,kCenters,kCentReScaled]=mlLabels_kmeans(dataForML,numLabel,nRow,nCol,lims,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

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