function [leftInds,rightInds,outRegions] = getSpecBounds(distBWorig,powerFilt,sampleNum,dupSpec)
% Get left and right spectrum boundaries

searchPoint=sampleNum*floor(dupSpec/2);
outRegions=zeros(size(distBWorig));

% Check if folding occurs
if length(find(~isnan(powerFilt(:,searchPoint))))<3
    leftInds=repmat(searchPoint,size(distBWorig,1),1);
    rightInds=leftInds+sampleNum-1;
    return
end

% Make sure regions don't touch
redoRegions=1;
distBWorig=double(distBWorig);
distBW=distBWorig;

while redoRegions

    regions=bwconncomp(distBW);

    redoRegions=0;
    for ii=1:regions.NumObjects
        regMask=zeros(size(distBW));
        regMask(regions.PixelIdxList{ii})=1;
        sumReg=any(regMask==1,1);
        sumSum=sum(sumReg);
        if sumSum>size(distBW,2)*0.8 % Connected
            distBW(regions.PixelIdxList{ii})=0;
            redoRegions=1;
            oneReg=1;
            while oneReg==1
                regMask=imerode(regMask,strel('disk',1));
                regions1=bwconncomp(regMask);
                oneReg=regions1.NumObjects;
            end
            distBW(regMask==1)=1;
        end
    end
end

distBWfirst=distBW;
distBWfirst(:,1:searchPoint-round(sampleNum/2))=nan;
distBWfirst(:,searchPoint+round(sampleNum/2):end)=nan;

if sum(distBWfirst(1,:),'omitnan')==0
    findFirstRow=min(find(sum(distBWfirst,2,'omitnan')>0));
    distBW(1,:)=distBWfirst(findFirstRow,:);
    regions=bwconncomp(distBW);
end

% Create gaps between clouds
powerSum=sum(powerFilt,2,'omitnan');
for jj=2:length(powerSum)
    if powerSum(jj)==0
        distBW(jj,:)=0;
    end
end

regions=bwconncomp(distBW);

% Find region closest to left boundary of middle spectrum

firstLine=distBW(1,:);
goodInds=find(firstLine==1);
if isempty(goodInds)
    error('No regions in first line.')
end

distToPoint=abs(goodInds-searchPoint);

[minNotUsed,minInd]=min(distToPoint);

regionPoint=sub2ind(size(distBW),1,goodInds(minInd));

ii=1;
regNum=[];

while isempty(regNum)
    regInds=regions.PixelIdxList{ii};
    inregion=intersect(regInds,regionPoint);

    if ~isempty(inregion)
        regNum=ii;
        outRegions(regInds)=1;
    end
    ii=ii+1;
end

sumCols=sum(outRegions,2);
firstEmpty=min(find(sumCols==0));

% Repeat procedure until we reach the last gate

while firstEmpty<=size(distBW,1)
    prevLine=outRegions(firstEmpty-1,:);
    searchPointL=min(find(prevLine==1));
    searchPointR=max(find(prevLine==1));

    jj=0;
    goodInds=[];
    checkLine=-1;
    while isempty(goodInds) & checkLine<size(distBW,1)
        checkLine=firstEmpty+jj;
        firstLine=distBW(checkLine,:);
        if jj==15
            searchPointL=searchPoint;
            searchPointR=searchPoint;
        end
        firstLine(1:searchPointL-sampleNum)=0;
        firstLine(searchPointR+sampleNum:end)=0;
        goodInds=find(firstLine==1);
        sumCols(checkLine)=1;
        jj=jj+1;
    end

    if isempty(goodInds)
        firstEmpty=size(distBW,1)+1;
        continue
    end

    distToPointL=abs(goodInds-searchPointL);
    [minL,minIndL]=min(distToPointL);
    distToPointR=abs(goodInds-searchPointR);
    [minR,minIndR]=min(distToPointR);

    if minL>=minR
        regionPoint=sub2ind(size(distBW),checkLine,goodInds(minIndR));
    else
        regionPoint=sub2ind(size(distBW),checkLine,goodInds(minIndL));
    end

    ii=1;
    regNum=[];

    while isempty(regNum)
        regInds=regions.PixelIdxList{ii};
        inregion=intersect(regInds,regionPoint);

        if ~isempty(inregion)
            regNum=ii;
            outRegions(regInds)=1;
        end
        ii=ii+1;
    end

    sumCols=sumCols+sum(outRegions,2);
    firstEmpty=min(find(sumCols==0));
end

% Find boundaries

leftInds=nan(size(distBW,1),1);

for ii=1:size(distBW,1)
    oneInds=find(outRegions(ii,:)==1);
    if isempty(oneInds)
        leftInds(ii)=leftInds(ii-1);
    else
        leftInds(ii)=round(mean(oneInds));
    end
end

rightInds=leftInds+sampleNum-1;

end