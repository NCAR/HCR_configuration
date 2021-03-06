function waterMasked = joinCloudParts(waterShed,refl)
% Joind areas that have a long border
ridgeLines=waterShed==0;
%ridgeLines(refl==0)=0;

ridgeSize=bwconncomp(ridgeLines);
ridges=ridgeSize.PixelIdxList;
waterMasked=waterShed;
waterMasked(refl==0)=0;

for kk=1:length(ridges)
    thisRidge=ridges{kk};
    if length(thisRidge)>50
        waterMasked(ridges{kk})=1;
        continue
    end
    thisRidgeIm=zeros(size(refl));
    thisRidgeIm(thisRidge)=1;
    largerRidge=imdilate(thisRidgeIm, strel('disk', 1));
    areaPix=waterMasked(largerRidge==1);
    areaPix(areaPix==0)=[];
    unPix=unique(areaPix);
    
    if length(unPix)==1
        waterMasked(ridges{kk})=1;
        continue
    end
    
    if length(unPix)>2
        disp(num2str(length(unPix)));
    end
    
    % Check circumfirence and proportions of adjacent regions
    cir=[];
    areaWidth=[];
    areaHeight=[];
    areaSize=[];
    for ll=1:length(unPix)
        % Mask surrounding areas
        maskUn=zeros(size(refl));
        maskUn(waterMasked==unPix(ll))=1;
        areaSize=[areaSize length(find(maskUn==1))];
        
        % Circumference
        cirAll=bwconncomp(maskUn);
        cir=[cir length(cirAll.PixelIdxList{1})];
        
        % Width
        widthAll=zeros(size(maskUn,1),1);
        for ii=1:size(maskUn,1)
            oneData=find(maskUn(ii,:)==1);
            if ~isempty(oneData)
                widthAll(ii)=max(oneData)-min(oneData);
            end
        end
        areaWidth=[areaWidth max(widthAll)];
        
        % Height
        heightAll=zeros(size(maskUn,2),1);
        for ii=1:size(maskUn,2)
            oneData=find(maskUn(:,ii)==1);
            if ~isempty(oneData)
                heightAll(ii)=max(oneData)-min(oneData);
            end
        end
        areaHeight=[areaHeight max(heightAll)];
    end
    maxFrac=max(length(ridges{kk})./cir);
    maxWidth=max(areaWidth);
    maxHeight=max(areaHeight);
    minSize=min(areaSize);
    if minSize<=10000 | maxFrac>0.001 | length(ridges{kk})/maxWidth>0.3 | length(ridges{kk})/maxHeight>0.3
        waterMasked(ridges{kk})=1;
    end
end
end

