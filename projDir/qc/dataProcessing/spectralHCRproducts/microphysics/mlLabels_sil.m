function [labelsInNum,kCenters,kCentReScaled]=mlLabels_sil(dataForML,numLabel,nRow,nCol,lims,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% k-means labels
disp('k-means');
%[kInds,kCenters]=kmeans(dataGood,numLabel,'MaxIter',1000);

array_silhoutte = zeros(1,num_kmeans);
for jj=1:numLabel
    % [cluster_assignments,centroids] = kmeans(X,j,'Distance','sqeuclidean','Start','sample');
    [kInds,kCenters]=kmeans(dataGood,numLabel,'MaxIter',1000);
    avgDWithin=zeros(dataset_len,1);
    avgDBetween=Inf(dataset_len,jj);
    for i=1:dataset_len
        for jj=1:jj
            boo=kInds==kInds(i);
            Xsamecluster=X(boo,:);
            if size(Xsamecluster,1)>1
                avgDWithin(i)=sum(sum((X(i,:)-Xsamecluster).^2,2))/(size(Xsamecluster,1)-1);
            end
            boo1= kInds~=kInds(i);
            Xdifferentcluster=X(boo1 & kInds ==jj,:);
            if ~isempty(Xdifferentcluster)
                avgDBetween(i,jj)=mean(sum((X(i,:)-Xdifferentcluster).^2,2));
            end
        end
    end
    % Calculate the silhouette values
    minavgDBetween = min(avgDBetween, [], 2);
    silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
    array_silhoutte(jj)=mean(silh);
    disp(num2str(array_silhoutte(jj)));
end
disp("Criterion values computed manually")
disp(array_silhoutte)

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

% Re-scale centers
kCentReScaled=nan(size(kCenters));
for ii=1:length(vars)
    kCentReScaled(:,ii)=kCenters(:,ii).*(lims.(vars{ii})(2)-lims.(vars{ii})(1))+lims.(vars{ii})(1);
end
end