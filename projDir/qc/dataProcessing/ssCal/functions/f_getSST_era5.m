% Reads era5 sea surface temperature data
function [SST]= f_getSST_era5(indir,inlon,inlat,intime)

inlon=wrapTo360(inlon);
refTime=datetime(1900,1,1,0,0,0);

%Initialize output
SST=nan(size(inlon));

%% Surface data
% Find ecmwf files
roundTime=dateshift(intime(end), 'start', 'hour', 'nearest');
monthStr=datestr(roundTime,'yyyymm');

sstFiles=dir([indir,'*.SSTK.*',monthStr,'0100_*.nc']);

if size(sstFiles,1)==0
    return
end

%read in lat and lon data
lonRean=ncread([sstFiles(1).folder,'/',sstFiles(1).name],'longitude');
latRean=ncread([sstFiles(1).folder,'/',sstFiles(1).name],'latitude');

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

for ii=1:size(ulonlatInd,1)
    timeReanS=ncread([sstFiles.folder,'/',sstFiles.name],'time');
    timeActualS=refTime+hours(timeReanS);
    timeIndS=find(timeActualS==roundTime);
    SST(ic==ii)=ncread([sstFiles.folder,'/',sstFiles.name],'SSTK',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1])-273.15;
end
end

