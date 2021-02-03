function stratConvSub=f_embeddedIsolated(stratConv)
% Subcategories are defined for each convective region by looking at
% the border
% 10: isolated tethered convective -> shares less than 10% of border with stratiform
% 11: embedded tethered convective-> shares more than 10% border with stratiform
% 12: isolated elevated convective -> shares less than 10% of border with stratiform
% 13: embedded elevated convective-> shares more than 10% border with stratiform

% For stratiform we look at the whole cloud
% 20: all stratiform or no embedded close by
% 21: stratiform with embedded convection

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

% % Embedded tethered convective
% if any(stratConv(cloudInds)==1 & any(stratConv(cloudInds)==2
%     stratConvSub(cloudInds)=10;
% end
% % Embedded elevated convective
% if any(stratConv(cloudInds)==0 & any(stratConv(cloudInds)==2
%     stratConvSub(cloudInds)=10;
%     return
% end

% All stratiform
if all(stratConv(cloudInds)==2)
    stratConvSub(cloudInds)=20;
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

convLarge=imdilate(convMask, strel('disk', 300));
convLarge= ~bwareaopen(~convLarge,50000);

% Strat with embedded
stratConvSub(stratConv==2 & convLarge==1)=21;
% Strat only
stratConvSub(stratConv==2 & convLarge==0)=20;
end