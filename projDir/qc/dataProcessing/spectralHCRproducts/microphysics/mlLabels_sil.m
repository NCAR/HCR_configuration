function [labelsInNum,kCenters,kCentReScaled]=mlLabels_sil(dataForML,numLabel,nRow,nCol,lims,vars)

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% k-means labels
dataset_len = size(dataGood,1);
array_silhoutte = zeros(1,numLabel);
for kk=2:numLabel
    disp(['k-means ',num2str(kk)]);
    disp(datetime)
    % [cluster_assignments,centroids] = kmeans(X,j,'Distance','sqeuclidean','Start','sample');
    [kInds,kCenters]=kmeans(dataGood,kk,'MaxIter',1000,'Replicates',5,'Display','final');
    avgDWithin=zeros(dataset_len,1);
    avgDBetween=Inf(dataset_len,kk);
    disp(['Silhouette ',num2str(kk)]);
    disp(datetime)
  
    % Get indices
    boos=nan(size(dataGood,1),kk);
    for jj=1:kk
        boos(:,jj)=kInds==jj;
    end
    for jj=1:kk
        % Same cluster
        boo=boos(:,jj);
        Xsamecluster=dataGood(boo==1,:);
        if size(Xsamecluster,1)>1
            avgDWithinThis=zeros(size(Xsamecluster,1),1);
            otherClus=1:kk;
            otherClus(jj)=[];
            avgDBetweenThis=inf(size(Xsamecluster,1),kk);
            for ii=1:size(Xsamecluster,1)
                avgDWithinThis(ii)=sum(sum((Xsamecluster(ii,:)-Xsamecluster).^2,2))/(size(Xsamecluster,1)-1);

                for ll=1:length(otherClus)
                    boo1=boos(:,otherClus(ll));
                    Xdifferentcluster=dataGood(boo1==1,:);
                    if ~isempty(Xdifferentcluster)
                        avgDBetweenThis(ii,otherClus(ll))=mean(sum((Xsamecluster(ii,:)-Xdifferentcluster).^2,2));
                    end
                end
                avgDWithin(boo==1,:)=avgDWithinThis;
                avgDBetween(boo==1,:)=avgDBetweenThis;
            end
        end
    end

    % Calculate the silhouette values
    minavgDBetween = min(avgDBetween, [], 2);
    silh = (minavgDBetween - avgDWithin) ./ max(avgDWithin,minavgDBetween);
    array_silhoutte(kk)=mean(silh);
    disp(['Value for ',num2str(kk),': ',num2str(array_silhoutte(kk))]);
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