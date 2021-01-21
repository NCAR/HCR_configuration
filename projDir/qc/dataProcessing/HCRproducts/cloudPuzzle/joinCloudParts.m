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

%% Loop through areas and if they are small, join at the smallest ridge
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
    
    if length(pixList)==0
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
            
            largerRidges(labelRidges==sortBL(1,1))=1;
            waterMasked(labelRidges==sortBL(1,1))=1;
            usedBorders(labelRidges==sortBL(1,1))=1;
        end
    end
end
%% Large ridges
forLabel=bwconncomp(waterMasked);
maskLabel=labelmatrix(forLabel);

ridgesAll=bwconncomp(largerRidges);
ridges=ridgesAll.PixelIdxList;

for kk=1:ridgesAll.NumObjects
    thisRidge=ridges{kk};
    
    % Find bordering regions
    maskThis=zeros(size(largerRidges));
    maskThis(thisRidge)=1;
    
    % Skeleton
    thisRidgeSkel=bwskel(logical(maskThis));
    ridgePixSum=sum(sum(thisRidgeSkel));
    
    if ridgePixSum>100
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
        continue
    end
    
    surrPix=maskLabel(thisRidge);
    surrPix(surrPix==0)=[];
    unPix=unique(surrPix);
    
    if length(unPix)==1
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
        continue
    end
    
    % Check circumfirence
    cir=[];
    
    for ll=1:length(unPix)
        % Mask surrounding areas
        maskUn=zeros(size(waterShed));
        maskUn(maskLabel==unPix(ll))=1;
                
        % Circumference
        %cirAll=bwconncomp(maskUn);
        %[r,c]=ind2sub(size(largerRidges),find(maskUn==1));
        %cirThis=boundary([c,r]);
        %cirThis=convhull(c,r);
        cirThis=regionprops(maskUn,'Perimeter');
        cir=[cir cirThis.Perimeter];
    end
    
    if length(unPix)>2
        error('Ridge divides more than one area.');
        disp(num2str(length(unPix)));
    end
    
    maxFrac=max(ridgePixSum./cir);
    if maxFrac>0.02
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
    end
end
 
%% In the second round, we only take care of areas that are too small
% ridgesAll=bwconncomp(largerRidges);
% ridges=ridgesAll.PixelIdxList;
% 
% % Label areas with individual numbers
% forLabel=bwconncomp(waterMasked);
% maskLabel=labelmatrix(forLabel);            
% 
% for kk=1:ridgesAll.NumObjects
%     thisRidge=ridges{kk};
%             
%     % Find bordering regions
%     maskThis=zeros(size(largerRidges));
%     maskThis(thisRidge)=1;
%     
%     % Skeleton
%     thisRidgeSkel=bwskel(logical(maskThis));    
%     ridgePixSum=sum(sum(thisRidgeSkel));
%     
%     if ridgePixSum>100
%         waterMasked(ridges{kk})=1;
%         continue
%     end
%         
%     surrPix=maskLabel(thisRidge);
%     surrPix(surrPix==0)=[];
%     unPix=unique(surrPix);    
%        
%     if length(unPix)==1
%         waterMasked(ridges{kk})=1;
%         continue
%     end
%     
%     % Check size of adjacent regions
%     areaSize=[];
%     
%     for ll=1:length(unPix)
%         % Mask surrounding areas
%         maskUn=zeros(size(waterShed));
%         maskUn(maskLabel==unPix(ll))=1;
%         areaSize=[areaSize length(find(maskUn==1))];
%     end
%     
%     if length(unPix)>2
%         error('Ridge divides more than one area.');
%         disp(num2str(length(unPix)));
%     end
%     
%     minSize=min(areaSize);
%     if minSize<=10000 % 10000 in original
%         waterMasked(ridges{kk})=1;
%     end
% end

end

