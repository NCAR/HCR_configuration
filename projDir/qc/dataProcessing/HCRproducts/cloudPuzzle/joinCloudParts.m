function waterMasked = joinCloudParts(waterShed)
% Join areas that have a long border or where one area is small, etc.

% Output
waterMasked=zeros(size(waterShed));
waterMasked(waterShed>0)=1;
waterMasked=bwareaopen(waterMasked,1000);

waterRidges=zeros(size(waterShed));
waterRidges(waterShed==0)=1;

% Clean it up
largerRidges=imclose(waterRidges, strel('disk', 7));
largerRidges=imdilate(largerRidges,strel('disk',5));

%largerRidgesKeep=largerRidges;

% % Make sure each ridge divides only two regions
% ridgesIn=bwconncomp(largerRidges);
% forLabel=bwconncomp(waterMasked);
% maskLabel=labelmatrix(forLabel);
% 
% for ii=1:ridgesIn.NumObjects
%     areasRidge=maskLabel(ridgesIn.PixelIdxList{ii});
%     [aCount,~] = groupcounts(areasRidge);
%     aCount(1)=[];
%     % If three areas
%     if length(aCount)>2
%         % Shrink down to skeleton
%         ridgeMask=zeros(size(maskLabel));
%         ridgeMask(ridgesIn.PixelIdxList{ii})=1;
%         ridgeSkel=bwskel(logical(ridgeMask));
%         
%         % Move along skeleton and remove all data where three areas are
%         % encountered
%         testLabel=nan(size(maskLabel));
%         testLabel(ridgeMask==1)=maskLabel(ridgeMask==1);
%         
%         [skelR,skelC]=find(ridgeSkel==1);
%         boxHalfSize=6;
%         for jj=1:length(skelR)
%             boxLabel=testLabel(skelR(jj)-boxHalfSize:skelR(jj)+boxHalfSize,skelC(jj)-boxHalfSize:skelC(jj)+boxHalfSize);
%             uBox=unique(boxLabel);
%             uBox(uBox==0)=[];
%             uBox(isnan(uBox))=[];
%             if length(uBox)>2
%                 waterMasked(skelR(jj)-boxHalfSize:skelR(jj)+boxHalfSize,skelC(jj)-boxHalfSize:skelC(jj)+boxHalfSize)=0;
%                 largerRidges(skelR(jj)-boxHalfSize:skelR(jj)+boxHalfSize,skelC(jj)-boxHalfSize:skelC(jj)+boxHalfSize)=0;
%             end
%         end
%     end
% end

%% Loop through areas and if they are small, join at the biggest ridge
someSmall=1;

while someSmall
    ridgesAll=bwconncomp(largerRidges);
    
    % Label Ridges
    labelRidges=labelmatrix(ridgesAll);
    
    % Label areas with individual numbers
    forLabel=bwconncomp(waterMasked);
    maskLabel=labelmatrix(forLabel);
    
    % Go from smallest to largest
    pixList=forLabel.PixelIdxList;
    [~,I] = sort(cellfun(@length,pixList));
    pixList = pixList(I);
    
    sizePix=cellfun(@length,pixList);
    pixList(sizePix>30000)=[];
    
    if length(pixList)<=1
        someSmall=0;
        continue
    end
    
    usedBorders=zeros(size(waterMasked));
    
    for ii=1:length(pixList)
        thisAreaPix=pixList{ii};
        
        % Find borders
        bordersA=labelRidges(thisAreaPix);
        uBorders=unique(bordersA);
        uBorders(uBorders==0)=[];
        
        if isempty(uBorders)
            continue
        end
        
        % Check if any of the borders has already been used
        testBorders=zeros(size(waterMasked));
        
        % Loop through ridges
        ridgeBL=[];
        for jj=1:length(uBorders)
            % Find bordering regions
            maskThis=zeros(size(largerRidges));
            maskThis(labelRidges==uBorders(jj))=1;
            testBorders(labelRidges==uBorders(jj))=usedBorders(labelRidges==uBorders(jj));
            
            % Skeleton
            thisRidgeSkel=bwskel(logical(maskThis));
            ridgePixSum=sum(sum(thisRidgeSkel));
            
            %borderLength=length(find(labelRidges==uBorders(jj)));
            ridgeBL=cat(1,ridgeBL,[uBorders(jj),ridgePixSum]);
        end
        
        if sum(sum(testBorders))==0 % If no border has been used
            sortBL=sortrows(ridgeBL,2,'descend');
            
            largerRidges(labelRidges==sortBL(1,1))=0;
            waterMasked(labelRidges==sortBL(1,1))=1;
            usedBorders(labelRidges==sortBL(1,1))=1;
        end
    end
end

%% Loop through all areas and check if they share a very large border

ridgesAll=bwconncomp(largerRidges);

% Label Ridges
labelRidges=labelmatrix(ridgesAll);

% Label areas with individual numbers
forLabel=bwconncomp(waterMasked);
maskLabel=labelmatrix(forLabel);

% Go from smallest to largest
pixList=forLabel.PixelIdxList;

if length(pixList)==1
    return
end

for ii=1:length(pixList)
    thisAreaPix=pixList{ii};
    
    % Find borders
    bordersA=labelRidges(thisAreaPix);
    uBorders=unique(bordersA);
    uBorders(uBorders==0)=[];
    
    if isempty(uBorders)
        continue
    end
    
    % Loop through borders and keep length
    ridgeBL=[];
    for jj=1:length(uBorders)
        % Find bordering regions
        maskThis=zeros(size(largerRidges));
        maskThis(labelRidges==uBorders(jj))=1;
        
        % Skeleton
        thisRidgeSkel=bwskel(logical(maskThis));
        ridgePixSum=sum(sum(thisRidgeSkel));
        
        ridgeBL=cat(1,ridgeBL,[uBorders(jj),ridgePixSum]);
    end
       
    % Circumference of area
    % Mask surrounding areas
    maskUn=zeros(size(waterShed));
    maskUn(thisAreaPix)=1;
    
    % Circumference
    % Close to smooth boundary
    maskUnSm=imclose(maskUn,strel('disk',10));
    cirThis=regionprops(maskUnSm,'Perimeter');
    
    % Fraction
    borderFrac=double(ridgeBL(:,2))./cirThis.Perimeter;
    [maxFrac,maxInd]=max(borderFrac);
    maxLength=ridgeBL(maxInd,2);
    
    if maxFrac>0.02 | maxLength>15
        waterMasked(labelRidges==ridgeBL(maxInd,1))=1;
    end
end

end

