function [leftInds rightInds, outRegions] = getSpecBounds(distBW,sampleNum,dupSpec,leftPrev)
% Get left and right spectrum boundaries

regions=bwconncomp(distBW);

redoRegions=0;

% Make sure regions don't touch
for ii=1:regions.NumObjects
    regMask=zeros(size(distBW));
    regMask(regions.PixelIdxList{ii})=1;
    sumReg=any(regMask==1,1);
    sumSum=sum(sumReg);
    if sumSum>size(distBW,2)*0.8; % Connected
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

if redoRegions
    regions=bwconncomp(distBW);
end

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
    checkLine=-1;
    while isempty(goodInds) & checkLine<size(distBW,1)
        checkLine=firstEmpty+jj;
        firstLine=distBW(checkLine,:);
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