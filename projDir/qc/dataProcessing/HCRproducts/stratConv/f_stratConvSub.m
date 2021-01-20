function stratConvSub=f_stratConvSub(stratConv,meltingLayer)
% Subcategories are defined for each convective region by looking at
% the border
% 10: isolated convective -> shares less than 10% of border with stratiform
% 11: warm embedded -> shares more than 10% border with stratiform and is
% warm
% 12: cold embedded -> as 11 but cold
% 13: warm small -> small and majority warm
% 14: cold small -> small and majority cold

% For stratiform we look at the whole cloud
% 20: all stratiform or no embedded
% 21: stratiform with embedded convection (11 or 12)
% 22: cold small -> small and majority cold
% 23: warm small -> small and majority warm

stratConvSub=nan(size(stratConv));

cloudInds=find(~isnan(stratConv));

% All convective
if sum(stratConv(cloudInds)==1)>length(cloudInds)*0.95
    stratConvSub(cloudInds)=10;
    return
end
% All stratiform
if sum(stratConv(cloudInds)==2)>length(cloudInds)*0.95
    stratConvSub(cloudInds)=20;
    return
end

% All stratiform are now stratiform embedded
stratConvSub(stratConv==2)=21;

% Divide convective into warm and cold embedded
stratConvSub(stratConv==1 & meltingLayer==20)=12;
stratConvSub(stratConv==1 & meltingLayer==10)=11;

end