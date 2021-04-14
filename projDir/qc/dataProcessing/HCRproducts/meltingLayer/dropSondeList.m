function [dropList dropTimes]=dropSondeList(startTime,endTime,indir,dropFormat)
% Make list of dropsonde files and dates
infiles=dir([indir,'*.',dropFormat]);

name1=infiles(1).name;
if strcmp(dropFormat,'eol')
    startDateInd=2;
else
    startDateInd=length(name1)-17;
end

datesIn=[];
for ii=1:length(infiles)
    fileName=infiles(ii).name;
    datesIn=cat(1,datesIn,datetime(str2num(fileName(startDateInd:startDateInd+3)),...
        str2num(fileName(startDateInd+4:startDateInd+5)),str2num(fileName(startDateInd+6:startDateInd+7)),...
        str2num(fileName(startDateInd+9:startDateInd+10)),str2num(fileName(startDateInd+11:startDateInd+12)),...
        str2num(fileName(startDateInd+13:startDateInd+14))));
end

rightInds=find(datesIn>=startTime & datesIn<=endTime);
dropList=infiles(rightInds);
dropTimes=datesIn(rightInds);
end

