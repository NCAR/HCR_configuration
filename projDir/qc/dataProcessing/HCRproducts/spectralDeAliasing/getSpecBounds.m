function [leftInds rightInds] = getSpecBounds(distBW,sampleNum,dupSpec,leftPrev)
% Get left and right spectrum boundaries

regions=bwconncomp(distBW);

outRegions=zeros(size(distBW));

% Find region closest to left boundary of middle spectrum
searchPoint=sampleNum*floor(dupSpec/2);

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
    while isempty(goodInds)
        checkLine=firstEmpty+jj;
        firstLine=distBW(checkLine,:);
        goodInds=find(firstLine==1);
        sumCols(checkLine)=1;
        jj=jj+1;
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

% Shrink regions and find boundaries

leftInds=nan(size(distBW,1),1);

for ii=1:size(distBW,1)
    oneInds=find(outRegions(ii,:)==1);
    if isempty(oneInds)
        leftInds(ii)=leftInds(ii-1);
    else
        leftInds(ii)=round(mean(oneInds));
    end
end

% plot(leftInds)
% hold on
% plot(leftPrev)
% hold off

rightInds=leftInds+sampleNum-1;

end