function stratConv=f_stratConvPart_elev(textPart,dbzPart,stratConvThresh,dbzThresh,asl,topo)
% Partition stratiform and convective bases on reflectivity texture field
stratConv=nan(size(textPart));

% Separate purely by threshold
convOnly=zeros(size(textPart));
stratOnly=zeros(size(textPart));

convOnly(textPart>=stratConvThresh)=1;
stratOnly(textPart<stratConvThresh)=1;
stratOnly(dbzPart<dbzThresh)=1;

% Enlarge convective areas so areas that are close together grow
% First make sure we are not enlarging narrow horizontal bands
convCleaned=zeros(size(textPart));
for ii=1:size(convOnly,2)
    convCol=convOnly(:,ii);
    convCleaned(:,ii)=bwareaopen(convCol,3);
end

convLarge=imdilate(convCleaned, strel('disk', 30));
stratOnly(convLarge==1)=0;
stratOnly(convOnly==1 & convLarge==0)=1; % Stratiform

% Remove stratiform areas that are small: less than a certain number
% of pixels and percentage of total cloud area
percMin=0.2;

stratAreaThresh=round(sum(sum(~isnan(dbzPart)))*percMin);
stratCleaned=bwareaopen(stratOnly,stratAreaThresh);

stratConv(stratCleaned==1)=2; % Stratiform
stratConv(~isnan(dbzPart) & stratCleaned~=1)=0; % Convective

stratConv(isnan(stratConv) & ~isnan(dbzPart))=2;
stratConv(isnan(dbzPart))=nan;

% Loop through convective areas and check if they are near the surface
convMask=zeros(size(stratConv));
convMask(stratConv==0)=1;

convAreas=bwconncomp(convMask);

% Calculate distance between asl and topo
distAslTopo=asl-topo;

for ii=1:convAreas.NumObjects
    aslArea=asl(convAreas.PixelIdxList{ii});
    nearSurfPix=sum(aslArea<500);
    if nearSurfPix~=0
        stratConv(convAreas.PixelIdxList{ii})=1;
    end
end

% If 95% of total cloud is convective, the whole cloud is convective
if (sum(sum(stratConv==1))+sum(sum(stratConv==0)))>sum(sum(~isnan(dbzPart)))*0.95
    if sum(sum(stratConv==1))>0
        stratConv(~isnan(dbzPart))=1;
    else
        stratConv(~isnan(dbzPart))=0;
    end
end

% Remove small convective blobs
pixMin=10000;
convAreaThresh=min([pixMin,round(sum(sum(~isnan(dbzPart)))*percMin)]);

convMask=zeros(size(textPart));
convMask(stratConv==1)=1;
convMask(stratConv==0)=1;
convCleaned=bwareaopen(convMask,convAreaThresh);

stratConv(convCleaned==0 & convMask==1)=2;
stratConv(isnan(dbzPart))=nan;
end

