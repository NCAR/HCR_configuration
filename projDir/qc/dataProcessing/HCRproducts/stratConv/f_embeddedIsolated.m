function stratConvSub=f_embeddedIsolated(stratConv)
% Subcategories are defined for each convective region by looking at
% the border
% 10: isolated tethered convective -> shares less than 10% of border with stratiform
% 11: embedded tethered convective-> shares more than 10% border with stratiform
% 12: isolated elevated convective -> shares less than 10% of border with stratiform
% 13: embedded elevated convective-> shares more than 10% border with stratiform

% For stratiform we look at the whole cloud
% 20: no embedded close by
% 21: stratiform with embedded convection
% 22: all stratiform (this is changed to 20 in the main code)

stratConvSub=nan(size(stratConv));

cloudInds=find(~isnan(stratConv));

% Isolated tethered convective
if all(stratConv(cloudInds)==1)
    stratConvSub(cloudInds)=10;
    return
end

% Isolated elevated convective
if all(stratConv(cloudInds)==0)
    stratConvSub(cloudInds)=12;
    return
end

% All stratiform
if all(stratConv(cloudInds)==2)
    stratConvSub(cloudInds)=22;
    return
end

% If some are stratifom, just separate tethered and elevated
% Tethered
stratConvSub(stratConv==1)=11;
% Elevated
stratConvSub(stratConv==0)=13;

% Now we make convective really big to separate stratiform with embedded
% from stratiform only
convMask=zeros(size(stratConv));
convMask(stratConv==0 | stratConv==1)=1;

horLarge=imdilate(convMask, strel('line', 400,0));

convLarge=imdilate(convMask, strel('disk', 200));
convLarge=~bwareaopen(~convLarge,50000);
convLarge=imclose(convLarge,strel('disk', 50));

convLarge(isnan(stratConv))=0;
convLarge=imfill(convLarge,'holes');

% Make sure we don't enlarge into unconnected areas
convRays=find(any(convLarge==1,1));
for ii=1:length(convRays)
    convCol=convLarge(:,convRays(ii));
    convHorCol=horLarge(:,convRays(ii));
    rayPieces=bwconncomp(convCol);
    for jj=1:rayPieces.NumObjects
        if ~any(convHorCol(rayPieces.PixelIdxList{jj})==1)
            convCol(rayPieces.PixelIdxList{jj})=0;
        end
    end
    convLarge(:,convRays(ii))=convCol;
end

% Strat with embedded
stratConvSub(stratConv==2 & convLarge==1)=21;
% Strat only
stratConvSub(stratConv==2 & convLarge==0)=20;

% Merge joining conv tethered and elevated
tethElev=zeros(size(stratConv));
tethElev(stratConvSub==11 | stratConvSub==13)=1;

tethElevAreas=bwconncomp(tethElev);
for ii=1:tethElevAreas.NumObjects;
    if any(stratConvSub(tethElevAreas.PixelIdxList{ii})==11)
        stratConvSub(tethElevAreas.PixelIdxList{ii})=11;
    end
end
end