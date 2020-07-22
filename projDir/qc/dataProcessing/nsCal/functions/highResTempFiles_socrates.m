%Find right file for high res temperature data
function [outFile]=highResTempFiles_socrates(startTime,endTime,indir)
allFileList=dir([indir,'*temperatures.txt']);

fileList={};
filenames={allFileList(:,:).name};
for ii=1:length(filenames)
    fileStart=datetime(str2num(filenames{ii}(15:18)),str2num(filenames{ii}(19:20)),...
        str2num(filenames{ii}(21:22)),str2num(filenames{ii}(24:25)),...
        str2num(filenames{ii}(26:27)),str2num(filenames{ii}(28:29)));
    fileEnd=datetime(str2num(filenames{ii}(34:37)),str2num(filenames{ii}(38:39)),...
        str2num(filenames{ii}(40:41)),str2num(filenames{ii}(43:44)),...
        str2num(filenames{ii}(45:46)),str2num(filenames{ii}(47:48)));
    if fileEnd>=startTime & fileStart<=endTime
        fileList{end+1}=[allFileList(ii).folder,'/',allFileList(ii).name];
    end
end

if length(fileList)==1
    outFile=fileList{1};
else
    disp('More than one temperature file found.');
end
end
