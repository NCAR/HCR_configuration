function dataLong=unJoinOverMissing(data,nonMissingInds)

dataInVars=fields(data);

dataLong=[];
for ii=1:length(dataInVars)
    if strcmp(dataInVars{ii},'time')
        stopHere=1;
        dataLong.time=data.time(nonMissingInds==1);
    else
        dataLong.(dataInVars{ii})=nan(size(data.(dataInVars{ii}),1),length(nonMissingInds));
        dataLong.(dataInVars{ii})(:,nonMissingInds==1)=data.(dataInVars{ii});
    end
end