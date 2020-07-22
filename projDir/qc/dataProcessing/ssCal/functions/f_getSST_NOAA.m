% Reads era5 sea surface temperature data
function [SST]= f_getSST_era5(indir,inlon,inlat,intime,ulonlatInd)
inyear=datestr(intime(1),'YYYY');
inday=datetime(year(intime(1)),month(intime(1)),day(intime(1)),0,0,0);

sstFile=dir([indir,'sst.day.mean.',inyear,'.nc']);

% inlon=wrapTo360(inlon);
refTime=datetime(1800,1,1,0,0,0);

%Initialize output
SST=nan;

%% Surface data
% Find file
if size(sstFile,1)==0
    return
end

% read in time
timeRean=ncread([sstFile(1).folder,'/',sstFile(1).name],'time');
timeActual=refTime+days(timeRean);

timeInd=find(timeActual==inday);

%read in lat and lon data
lonRean=ncread([sstFile(1).folder,'/',sstFile(1).name],'lon');
latRean=ncread([sstFile(1).folder,'/',sstFile(1).name],'lat');

lonInd=nan(size(inlon));
latInd=nan(size(inlat));
for ii=1:length(lonInd)
    lonIndAll=find(abs(lonRean-inlon(ii))==min(abs(lonRean-inlon(ii))));
    lonInd(ii)=lonIndAll(1);
    latIndAll=find(abs(latRean-inlat(ii))==min(abs(latRean-inlat(ii))));
    latInd(ii)=latIndAll(1);
end

lonlatInd=cat(2,lonInd,latInd);
[ulonlatInd,ia,ic]=unique(lonlatInd,'rows');

SSTall=[];
for ii=1:size(ulonlatInd,1)
    SSTall=[SSTall,ncread([sstFile.folder,'/',sstFile.name],'sst',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeInd],[1,1,1])];
end
SST=nanmean(SSTall);
end

