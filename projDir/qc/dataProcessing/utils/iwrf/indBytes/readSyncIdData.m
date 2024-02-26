function [syncIDthis,fileID]=readSyncIdData(fileID)

syncIDthis=fread(fileID,[1,2],'int32');

end