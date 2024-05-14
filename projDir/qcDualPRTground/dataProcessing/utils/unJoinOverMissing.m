function dataLong=unJoinOverMissing(data,nonMissingInds)

dataInVars=fields(data);

dataLong=[];
for ii=1:length(dataInVars)
    if strcmp(dataInVars{ii},'time')
        dataLong.time=NaT(size(data.(dataInVars{ii}),1),length(nonMissingInds));
        dataLong.time(:,nonMissingInds==1)=data.time;
    else
        dataLong.(dataInVars{ii})=nan(size(data.(dataInVars{ii}),1),length(nonMissingInds));
        dataLong.(dataInVars{ii})(:,nonMissingInds==1)=data.(dataInVars{ii});
    end
end