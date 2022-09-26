function [dataShort,nonMissingInds]=joinOverMissing(data,gapSecs)

nonMissingInds=findNonMissingInds(data,gapSecs);

dataInVars=fields(data);

dataShort=[];
for ii=1:length(dataInVars)
    dataShort.(dataInVars{ii})=data.(dataInVars{ii})(:,nonMissingInds==1);
end

end