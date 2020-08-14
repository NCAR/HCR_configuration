function cloudParams=calcCloudParams(dataCut);
% Calculate parameters for all individual clouds

cloudParams=[];

%% Max and min altitude and precipitation

% min/max asl
minAslAll=nan(1,size(dataCut.cloudPuzzle,2));
maxAslAll=nan(1,size(dataCut.cloudPuzzle,2));
precip=nan(1,size(dataCut.cloudPuzzle,2));

for jj=1:size(dataCut.cloudPuzzle,2)
    % Min/max alt
    aslCol=dataCut.asl(:,jj);
    [minAslAll(jj) minIndAsl]=min(aslCol);
    [maxAslAll(jj) maxIndAsl]=max(aslCol);
    % Check if flying in cloud
    % Pointing down
    if dataCut.elevation(jj)<0 & maxAslAll(jj)<10000 & maxIndAsl==18
        maxAslAll(jj)=nan;
    end
    % Pointing up
    if dataCut.elevation(jj)>=0 & minAslAll(jj)<10000 & minIndAsl==18
        minAslAll(jj)=nan;
    end
    
    % Precip
    if dataCut.elevation(jj)<0
        surfInd=min(find(dataCut.FLAG(:,jj)==7));
        if ~isempty(surfInd)
            precip(jj)=mean(dataCut.cloudPuzzle(surfInd-5:surfInd-1,jj));
        elseif any(dataCut.flagCensored(:,jj)==3)
            precip(jj)=1;
        end
    end
end
percWanted=0.1;

% Make sure we have enough data and calculate percentiles
if length(find(~isnan(minAslAll)))>length(minAslAll)/2
    sortedMin=sort(minAslAll,'ascend');
    percIndMin=round(percWanted*length(minAslAll));
    minAsl=sortedMin(percIndMin);
else
    minAsl=nan;
end

if length(find(~isnan(maxAslAll)))>length(maxAslAll)/2
    sortedMax=sort(maxAslAll,'descend');
    percIndMax=round(percWanted*length(maxAslAll));
    maxAsl=sortedMax(percIndMax);
else
    maxAsl=nan;
end

% Above ground level
cloudParams.minAgl=minAsl-mean(dataCut.TOPO,'omitnan')./1000;
cloudParams.maxAgl=maxAsl-mean(dataCut.TOPO,'omitnan')./1000;

% Precip
if length(find(~isnan(precip)))>length(precip)*0.03
    precCloud=1;
else
    precCloud=0;
end

cloudParams.precip=precCloud;
end

