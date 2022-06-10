function era5data = read_era5_levelFiles_uv1D(indir,startTime,endTime,SSTyes)
% Find the right time span for era5 data and read them in

refTime=datetime(1900,1,1,0,0,0);

%Initialize output
era5data=[];

% % Get intime hours
% inhoursAll=datetime(year(intime),month(intime),day(intime),hour(intime),0,0);
% inhoursTemp=unique(inhoursAll);
% inhoursLast=inhoursTemp(end)+hours(1);
% 
% ll=1;
% inhours=inhoursTemp(1);
% while inhours(end)~=inhoursLast
%     inhours=[inhours inhoursTemp(1)+hours(ll)];
%     ll=ll+1;
% end

startHour=datetime(year(startTime),month(startTime),day(startTime),hour(startTime),0,0);
endHour=datetime(year(endTime),month(endTime),day(endTime),hour(endTime)+1,0,0);

inhours=startHour:hours(1):endHour;

pHours=[];
tHours=[];
rhHours=[];
zHours=[];
uHours=[];
vHours=[];
pSurf=[];
tSurf=[];
rhSurf=[];
zSurf=[];
uSurf=[];
vSurf=[];
if SSTyes
    sstHours=[];
    sstSurf=[];
end

for kk=1:length(inhours) %Loop through all hours
    %% Pressure level data
    % Find ecmwf files
    roundTime=inhours(kk);
    dayStr=datestr(roundTime,'yyyymmdd');
    
    rFiles=dir([indir,'R.*',dayStr,'00_',dayStr,'23.nc']);
    tFiles=dir([indir,'T.*',dayStr,'00_',dayStr,'23.nc']);
    zFiles=dir([indir,'Z.*',dayStr,'00_',dayStr,'23.nc']);
            
    if size(zFiles,1)==0 | size(zFiles,1)==0 | size(zFiles,1)==0
        disp('No model data found.');
        return
    end
    
    %read in time, lat and lon data
    timeRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'time');
    timeActual=refTime+hours(timeRean);
    timeInd=find(timeActual==roundTime);
    
    lonRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'longitude');
    latRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'latitude');
    
    t=nan(length(lonRean),length(latRean),length(rFiles));
    rh=nan(length(lonRean),length(latRean),length(rFiles));
    z=nan(length(lonRean),length(latRean),length(rFiles));
    p=nan(length(lonRean),length(latRean),length(rFiles));
    
    for jj=1:size(tFiles,1)
        p(:,:,jj)=fliplr(ncread([rFiles(jj).folder,'/',rFiles(jj).name],'level'));
        t(:,:,jj)=fliplr(squeeze(ncread([tFiles(jj).folder,'/',tFiles(jj).name],'T',[1,1,1,timeInd],[inf,inf,inf,1]))-273.15);
        rh(:,:,jj)=fliplr(squeeze(ncread([rFiles(jj).folder,'/',rFiles(jj).name],'R',[1,1,1,timeInd],[inf,inf,inf,1])));
        z(:,:,jj)=fliplr(squeeze(ncread([zFiles(jj).folder,'/',zFiles(jj).name],'Z',[1,1,1,timeInd],[inf,inf,inf,1])));
    end
    
    z=z./9.806;
    
    pHours=cat(4,pHours,p);
    tHours=cat(4,tHours,t);
    rhHours=cat(4,rhHours,rh);
    zHours=cat(4,zHours,z);
    
    %% Surface data
    % Find ecmwf files
    monthStr=datestr(roundTime,'yyyymm');
    
    dFiles=dir([indir,'VAR_2D.*',monthStr,'0100_*.nc']);
    tsFiles=dir([indir,'VAR_2T.*',monthStr,'0100_*.nc']);
    uFiles=dir([indir,'VAR_10U.*',monthStr,'0100_*.nc']);
    vFiles=dir([indir,'VAR_10V.*',monthStr,'0100_*.nc']);
    pFiles=dir([indir,'SP.*',monthStr,'0100_*.nc']);
    
    info=ncinfo([dFiles.folder,'/',dFiles.name]);
    
    timeReanS=ncread([dFiles.folder,'/',dFiles.name],'time');
    timeActualS=refTime+hours(timeReanS);
    timeIndS=find(timeActualS==roundTime);
    pS=fliplr(ncread([pFiles.folder,'/',pFiles.name],'SP',[1,1,timeIndS],[inf,inf,1])./100);
    tS=fliplr(ncread([tsFiles.folder,'/',tsFiles.name],'VAR_2T',[1,1,timeIndS],[inf,inf,1])-273.15);
    td=fliplr(ncread([dFiles.folder,'/',dFiles.name],'VAR_2D',[1,1,timeIndS],[inf,inf,1])-273.15);
    rhS=100*(exp((17.625*td)./(243.04+td))./exp((17.625*tS)./(243.04+tS)));
    sfcU=fliplr(ncread([uFiles.folder,'/',uFiles.name],'VAR_10U',[1,1,timeIndS],[inf,inf,1]));
    sfcV=fliplr(ncread([vFiles.folder,'/',vFiles.name],'VAR_10V',[1,1,timeIndS],[inf,inf,1]));
    
    pSurf=cat(3,pSurf,pS);
    tSurf=cat(3,tSurf,tS);
    rhSurf=cat(3,rhSurf,rhS);
    uSurf=cat(3,uSurf,sfcU);
    vSurf=cat(3,vSurf,sfcV);
    
    %% SST
    if SSTyes
        sstFiles=dir([indir,'SSTK.*',monthStr,'0100_*.nc']);
        sst=fliplr(ncread([sstFiles.folder,'/',sstFiles.name],'SSTK',[1,1,timeIndS],[inf,inf,1])-273.15);
        sstSurf=cat(3,sstSurf,sst);
    end
    
end
era5data.lat=latRean;
era5data.lon=lonRean;
era5data.Temperature=tHours;
era5data.rh=rhHours;
era5data.z=zHours;
era5data.p=pHours;
era5data.pSurf=pSurf;
era5data.tSurf=tSurf;
era5data.rhSurf=rhSurf;
era5data.uSurf=uSurf;
era5data.vSurf=vSurf;
era5data.time=inhours;
if SSTyes
    era5data.sstSurf=sstSurf;
end

end
