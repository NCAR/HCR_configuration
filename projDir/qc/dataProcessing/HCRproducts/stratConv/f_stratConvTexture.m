function stratConv=f_stratConvTexture(textPart,dbzPart,stratConvThresh,dbzThresh)
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
    convCleaned(:,ii)=bwareaopen(convCol,10);
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
stratConv(~isnan(dbzPart) & stratCleaned~=1)=1; % Convective

stratConv(isnan(dbzPart))=nan;

% Remove small convective blobs
pixMin=10000;
convAreaThresh=min([pixMin,round(sum(sum(~isnan(dbzPart)))*percMin)]);

convMask=zeros(size(textPart));
convMask(stratConv==1)=1;
convCleaned=bwareaopen(convMask,convAreaThresh);

stratConv(convCleaned==0 & convMask==1)=2;
stratConv(isnan(dbzPart))=nan;

% Handle small clouds
smallThresh=20000;
if round(sum(sum(~isnan(dbzPart))))<smallThresh
    % Assign classification that is more than 50%
    convPix=length(find(stratConv==1));
    stratPix=length(find(stratConv==2));
    if convPix>=stratPix
       stratConv(~isnan(dbzPart))=1;
    else
        stratConv(~isnan(dbzPart))=2;
    end
    
    % If more than x% of dbz data is >x dBZ -> convective
    reflCut=5;
    reflPerc=0.03;
    if length(find(dbzPart>reflCut))/round(sum(sum(~isnan(dbzPart))))>reflPerc
        stratConv(~isnan(dbzPart))=1;
    end
end
end

