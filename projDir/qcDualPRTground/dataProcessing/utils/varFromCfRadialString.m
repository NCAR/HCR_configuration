function [outValue] = varFromCfRadialString(infile,varName)
%Read variable that is stored in cfradial status_xml string
status_xml_in=ncread(infile,'status_xml');
status_xml_lines=(strsplit(status_xml_in','\n'))';
try
    temp_lines=strfind(status_xml_lines,varName);
    tempInds=find(cellfun(@isempty,temp_lines)==0);
    tempCellsLines=status_xml_lines(tempInds);
    tempCellsSplit = regexp(tempCellsLines, '[><]', 'split');
    rightCell=find(cellfun(@(c) any(strcmp(c, varName)),tempCellsSplit));
    outValue=str2num(tempCellsSplit{rightCell,1}{1,3});
catch
    disp(['Variable ' varName ' not found. Passing NAN.']);
    outValue=nan;
end
end

