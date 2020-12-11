function [dropAlt,dropT]=getDropData(dropList,fileFormat)
% Get dropsonde altitude and temperature from file list
dropAlt={};
dropT={};
for ii=1:length(dropList)
    disp(dropList(ii).name);
    fileName=[dropList(ii).folder,'/',dropList(ii).name];
    
    % The netcdf code is not tested yet!!!!!!
    if strcmp(fileFormat,'nc')
        dropAlt{end+1}=ncread(fileName,'alt');
        dropT{end+1}=ncread(fileName,'tdry');
    else
        allVars=readtable(fileName,'FileType','text');
        tempIn=allVars.Var7;
        tempIn(tempIn==-999)=nan;
        dropT{end+1}=tempIn;
        altIn=allVars.Var17;
        altIn(altIn==-999)=nan;
        dropAlt{end+1}=altIn;
    end
end
end

