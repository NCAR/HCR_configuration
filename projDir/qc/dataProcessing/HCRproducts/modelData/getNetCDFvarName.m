function varName=getNetCDFvarName(infile,searchString)
% Find variable name in netcdf file
infoIn=ncinfo(infile);
varName=[];

ii=1;
while isempty(varName) & ii<=length(infoIn.Variables)
    varNamesCheck=infoIn.Variables(ii).Name;
    foundStr=strfind(varNamesCheck,searchString);
    if ~isempty(foundStr)
        varName=varNamesCheck;
    end
    ii=ii+1;
end
end