function classBasic=f_classBasic(convOne,stratMixed,mixedConv)
% Find basic classification and merge suitable data
classBasic=nan(size(convOne));

%% Handle mixed
% Make mixed+conv mask
maskMixed=convOne>=stratMixed;

% Check if all stratiform
if max(max(maskMixed))==0;
    classBasic(~isnan(convOne))=1;
    return
end

% Remove areas that are small
maskMixed=bwareaopen(maskMixed,500);

% Set up check
horLarge=imdilate(maskMixed, strel('line', 100,0));

% Enlarge mixedConv
mixedLarge1=imdilate(maskMixed, strel('disk', 25)); %30
mixedLarge=imclose(mixedLarge1,strel('disk', 50));
mixedLarge(isnan(convOne))=0;
mixedLarge=imfill(mixedLarge,'holes');
mixedLarge=imerode(mixedLarge,strel('disk', 3));

% Make sure we don't enlarge into unconnected areas
mixedRays=find(any(mixedLarge==1,1));
for ii=1:length(mixedRays)
    mixedCol=mixedLarge(:,mixedRays(ii));
    mixedHorCol=horLarge(:,mixedRays(ii));
    rayPieces=bwconncomp(mixedCol);
    if rayPieces.NumObjects>1
        for jj=1:rayPieces.NumObjects
            if ~any(mixedHorCol(rayPieces.PixelIdxList{jj})==1)
                mixedCol(rayPieces.PixelIdxList{jj})=0;
            end
        end
        mixedLarge(:,mixedRays(ii))=mixedCol;
    end
end

mixedLarge=imdilate(mixedLarge,strel('disk', 3));

%% Handle conv
% Make conv mask
maskConv=convOne>=mixedConv;

% Set up check
horLarge2=imdilate(maskConv, strel('line', 100,0));

% Enlarge conv
convLarge1=imdilate(maskConv, strel('disk', 15)); %30
convLarge=imclose(convLarge1,strel('disk', 50));
convLarge(isnan(convOne))=0;
convLarge=imfill(convLarge,'holes');
convLarge=imerode(convLarge,strel('disk', 3));

% Make sure we don't enlarge into unconnected areas
convRays=find(any(convLarge==1,1));
for ii=1:length(convRays)
    convCol=convLarge(:,convRays(ii));
    convHorCol=horLarge2(:,convRays(ii));
    rayPieces2=bwconncomp(convCol);
    if rayPieces2.NumObjects>1
        for jj=1:rayPieces2.NumObjects
            if ~any(convHorCol(rayPieces2.PixelIdxList{jj})==1)
                convCol(rayPieces2.PixelIdxList{jj})=0;
            end
        end
        convLarge(:,convRays(ii))=convCol;
    end
end

convLarge=imdilate(convLarge,strel('disk', 3));

classBasic(convLarge)=3; % Convective
classBasic(isnan(classBasic) & mixedLarge)=2; % Mixed
classBasic(isnan(classBasic))=1; % Stratiform
classBasic(isnan(convOne))=nan;

% If 80% of total cloud is convective, the whole cloud is convective
if sum(sum(convLarge))>sum(sum(~isnan(convOne)))*0.8
    classBasic(~isnan(convOne))=3;
    return
end

% If 70% of total cloud is convective or mixed, stratiform areas become
% mixed
if sum(sum(mixedLarge))>sum(sum(~isnan(convOne)))*0.7
    classBasic(classBasic==1)=2;
end

% Unconnected stratiform areas that are small are set to mixed
stratOnly=classBasic==1;
percMin=0.2;

stratAreaThresh=round(sum(sum(~isnan(convOne)))*percMin);
stratCleaned=bwareaopen(stratOnly,stratAreaThresh);

classBasic(classBasic==1 & stratCleaned==0)=2; % Mixed

% Take care of mixed borders
maskMixedEnd=classBasic==2;
maskMixedEnd=imdilate(maskMixedEnd, strel('disk', 1));

mixedAreas=bwconncomp(maskMixedEnd);

% If it has no neighboring stratiform but neighboring convective and is
% small: set to convective
for ii=1:mixedAreas.NumObjects
    classFeat=classBasic(mixedAreas.PixelIdxList{ii});
    if ~any(classFeat==1) & any(classFeat==3) & length(mixedAreas.PixelIdxList{ii})<sum(sum(~isnan(convOne)))*0.2
        classBasic(mixedAreas.PixelIdxList{ii})=3;
    end
end

classBasic(isnan(convOne))=nan;
end

