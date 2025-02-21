function [labelsInNum,labelsOptEll,labNumEll]=mlLabels(dataForML,numLabel,nRow,nCol)

labelsInNum=nan;
labelsOptEll=nan;
labNumEll=nan;

goodInds=~any(isnan(dataForML),2);
dataGood=dataForML(goodInds==1,:);

% Labels with pre-defined number of clusters
kInds=kmeans(dataGood,numLabel,'MaxIter',1000);

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kInds;

labelsInNum=reshape(kIndsAll,nRow,nCol);

% Labels with optimum number of clusters using ellbow method
disp('Ellbow optimization ...');
[kIndsOpt,C,SUMD,labNumEll]=kmeans_opt(dataGood,15);

kIndsAll=nan(size(dataForML,1),1);
kIndsAll(goodInds==1)=kIndsOpt;

labelsOptEll=reshape(kIndsAll,nRow,nCol);
end