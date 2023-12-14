function data=trimFirstFile(data,endInd)
dataFields=fields(data);

checkLength=length(data.time);

for ii=1:length(dataFields)
    thisField=data.(dataFields{ii});
    if size(thisField,2)==checkLength
        data.(dataFields{ii})(:,1:endInd)=[];
    else
        data=rmfield(data,dataFields{ii});
    end
end
end