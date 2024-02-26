function [statusXml,fileID]=readStatusXmlData(fileID)
statusXml.xmlLength=fread(fileID,1,'int32',17*4);
statusXml.text=fread(fileID,[1,statusXml.xmlLength-6],'*char');
end