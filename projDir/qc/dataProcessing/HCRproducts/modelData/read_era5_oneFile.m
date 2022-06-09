function era5data = read_era5_oneFile(indir,startTime,endTime,SSTyes)
% Find the right time span for era5 data and read them in

%Initialize output
era5data=[];

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
uSurf=[];
vSurf=[];
if SSTyes
    sstSurf=[];
end

for kk=1:length(inhours) %Loop through all hours
    %% Pressure level data
    % Find ecmwf files
    roundTime=inhours(kk);
    dayStr=datestr(roundTime,'yyyymmdd');
    
    rFile=dir([indir,'*_r.*',dayStr,'00_',dayStr,'23.*.nc']);
    tFile=dir([indir,'*_t.*',dayStr,'00_',dayStr,'23.*.nc']);
    zFile=dir([indir,'*_z.*',dayStr,'00_',dayStr,'23.*.nc']);
    uFile=dir([indir,'*_u.*',dayStr,'00_',dayStr,'23.*.nc']);
    vFile=dir([indir,'*_v.*',dayStr,'00_',dayStr,'23.*.nc']);
            
    if size(zFile,1)==0 | size(zFile,1)==0 | size(zFile,1)==0
        disp('No model data found.');
        return
    end
    % Get variable names
    timeName=getNetCDFvarName([rFile.folder,'/',rFile.name],'time');
    longName=getNetCDFvarName([rFile.folder,'/',rFile.name],'lon');
    latName=getNetCDFvarName([rFile.folder,'/',rFile.name],'lat');
    pName=getNetCDFvarName([rFile.folder,'/',rFile.name],'level');
    if isempty(pName)
        pName=getNetCDFvarName([rFile.folder,'/',rFile.name],'lv');
    end
    tName=getNetCDFvarName([tFile.folder,'/',tFile.name],'T');
    rName=getNetCDFvarName([rFile.folder,'/',rFile.name],'R');
    zName=getNetCDFvarName([zFile.folder,'/',zFile.name],'Z');
    uName=getNetCDFvarName([uFile.folder,'/',uFile.name],'U');
    vName=getNetCDFvarName([vFile.folder,'/',vFile.name],'V');
    
    %read in time, lat and lon data
    timeRean=ncread([rFile.folder,'/',rFile.name],timeName);
    if strcmp(timeName,'time')
        refTime=datetime(1900,1,1,0,0,0);
        timeActual=refTime+hours(timeRean);
    elseif strcmp(timeName,'forecast_time0')
        refTime=datetime(year(roundTime),month(roundTime),day(roundTime));
        timeActual=refTime+hours(timeRean);
    else
        warning('Something is wrong with the era5 time.')
    end
    timeInd=find(timeActual==roundTime);
    
    lonRean=ncread([rFile.folder,'/',rFile.name],longName);
    latRean=ncread([rFile.folder,'/',rFile.name],latName);
    pRean=ncread([rFile.folder,'/',rFile.name],pName);

    %p(:,:,jj)=fliplr(ncread([rFile.folder,'/',rFile.name],pName));
    t=fliplr(squeeze(ncread([tFile.folder,'/',tFile.name],tName,[1,1,1,timeInd],[inf,inf,inf,1]))-273.15);
    rh=fliplr(squeeze(ncread([rFile.folder,'/',rFile.name],rName,[1,1,1,timeInd],[inf,inf,inf,1])));
    z=fliplr(squeeze(ncread([zFile.folder,'/',zFile.name],zName,[1,1,1,timeInd],[inf,inf,inf,1])));
    u=fliplr(squeeze(ncread([uFile.folder,'/',uFile.name],uName,[1,1,1,timeInd],[inf,inf,inf,1])));
    v=fliplr(squeeze(ncread([vFile.folder,'/',vFile.name],vName,[1,1,1,timeInd],[inf,inf,inf,1])));
    p=repmat(pRean,1,length(lonRean),length(latRean));
    p=double(permute(p,[2,3,1]));

    z=z./9.806;
    
    pHours=cat(4,pHours,p);
    tHours=cat(4,tHours,t);
    rhHours=cat(4,rhHours,rh);
    zHours=cat(4,zHours,z);
    uHours=cat(4,uHours,u);
    vHours=cat(4,vHours,v);
    
    %% Surface data
    % Find ecmwf files
    monthStr=datestr(roundTime,'yyyymm');
    
    dFiles=dir([indir,'*_2d.*',monthStr,'0100_*.nc']);
    tsFiles=dir([indir,'*_2t.*',monthStr,'0100_*.nc']);
    usFiles=dir([indir,'*_10u.*',monthStr,'0100_*.nc']);
    vsFiles=dir([indir,'*_10v.*',monthStr,'0100_*.nc']);
    pFiles=dir([indir,'*_sp.*',monthStr,'0100_*.nc']);
    
    % Get var names
    timeNameS=getNetCDFvarName([dFiles.folder,'/',dFiles.name],'time');
    spNameS=getNetCDFvarName([pFiles.folder,'/',pFiles.name],'SP');
    tNameS=getNetCDFvarName([tsFiles.folder,'/',tsFiles.name],'2T');
    dNameS=getNetCDFvarName([dFiles.folder,'/',dFiles.name],'2D');
    uNameS=getNetCDFvarName([usFiles.folder,'/',usFiles.name],'10U');
    vNameS=getNetCDFvarName([vsFiles.folder,'/',vsFiles.name],'10V');

    % Get data
    timeReanS=ncread([dFiles.folder,'/',dFiles.name],timeNameS);
    timeActualS=refTime+hours(timeReanS);
    timeIndS=find(timeActualS==roundTime);
    pS=fliplr(ncread([pFiles.folder,'/',pFiles.name],spNameS,[1,1,timeIndS],[inf,inf,1])./100);
    tS=fliplr(ncread([tsFiles.folder,'/',tsFiles.name],tNameS,[1,1,timeIndS],[inf,inf,1])-273.15);
    td=fliplr(ncread([dFiles.folder,'/',dFiles.name],dNameS,[1,1,timeIndS],[inf,inf,1])-273.15);
    rhS=100*(exp((17.625*td)./(243.04+td))./exp((17.625*tS)./(243.04+tS)));
    sfcU=fliplr(ncread([usFiles.folder,'/',usFiles.name],uNameS,[1,1,timeIndS],[inf,inf,1]));
    sfcV=fliplr(ncread([vsFiles.folder,'/',vsFiles.name],vNameS,[1,1,timeIndS],[inf,inf,1]));
    
    pSurf=cat(3,pSurf,pS);
    tSurf=cat(3,tSurf,tS);
    rhSurf=cat(3,rhSurf,rhS);
    uSurf=cat(3,uSurf,sfcU);
    vSurf=cat(3,vSurf,sfcV);
    
    %% SST
    if SSTyes
        sstFiles=dir([indir,'*_sstk.*',monthStr,'0100_*.nc']);
        sstNameS=getNetCDFvarName([sstFiles.folder,'/',sstFiles.name],'SST');
        sst=fliplr(ncread([sstFiles.folder,'/',sstFiles.name],sstNameS,[1,1,timeIndS],[inf,inf,1])-273.15);
        sstSurf=cat(3,sstSurf,sst);
    end
    
end
era5data.lat=latRean;
era5data.lon=lonRean;
era5data.Temperature=tHours;
era5data.rh=rhHours;
era5data.z=zHours;
era5data.p=pHours;
era5data.u=uHours;
era5data.v=vHours;
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

