function dataThis=trimData(data,startInd,endInd)
dataFields=fields(data);

for ii=1:length(dataFields)
    thisField=data.(dataFields{ii});
    if size(thisField,2)==length(data.time)
        dataThis.(dataFields{ii})=thisField(:,startInd:endInd);
    else
        dataThis.(dataFields{ii})=thisField;
    end
end
end