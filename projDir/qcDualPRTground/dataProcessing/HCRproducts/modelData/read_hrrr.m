function hrrrData = read_hrrr(indir,startTime,endTime,SSTyes)
% Find the right time span for era5 data and read them in

refTime=datetime(1900,1,1,0,0,0);

%Initialize output
hrrrData=[];

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

lonlatFile=dir([indir,'netcdf/*.nc']);
lon=ncread([lonlatFile.folder,'/',lonlatFile.name],'gridlon_0');
hrrrData.lon=flipud(lon');
lat=ncread([lonlatFile.folder,'/',lonlatFile.name],'gridlat_0');
hrrrData.lat=flipud(lat');

for kk=1:length(inhours) %Loop through all hours
    %% Pressure level data
    % Find ecmwf files
    roundTime=inhours(kk);
    hourStr=datestr(roundTime,'HH');
    dayStr=datestr(roundTime,'yyyymmdd');
    
    hourFile=dir([indir,'/',dayStr,'/hrrr.t',hourStr,'z.wrfprsf00.grib2']);
    hourFile=[hourFile.folder,'/',hourFile.name];
               
    info=georasterinfo(hourFile);
    elem=info.Metadata.Element;
    descr=info.Metadata.Description;

    pressLevs=find(contains(descr,'ISBL'));
    m0levs=find(contains(descr,'0[-] SFC'));
    m2levs=find(contains(descr,'2[m] HTGL'));
    m10levs=find(contains(descr,'10[m] HTGL'));
    
    allT=find(elem=='TMP');
    bandTpress=ismember(allT,pressLevs);
    bandTpress=allT(bandTpress==1);
    bandTsurf=ismember(allT,m2levs);
    bandTsurf=allT(bandTsurf==1);

    allRH=find(elem=='RH');
    bandRHpress=ismember(allRH,pressLevs);
    bandRHpress=allRH(bandRHpress==1);
    bandRHsurf=ismember(allRH,m2levs);
    bandRHsurf=allRH(bandRHsurf==1);

    allZ=find(elem=='HGT');
    bandZpress=ismember(allZ,pressLevs);
    bandZpress=allZ(bandZpress==1);
    % bandZsurf=ismember(allZ,m2levs);
    % bandZsurf=allZ(bandZsurf==1);

    allU=find(elem=='UGRD');
    bandUpress=ismember(allU,pressLevs);
    bandUpress=allU(bandUpress==1);
    bandUsurf=ismember(allU,m10levs);
    bandUsurf=allU(bandUsurf==1);

    allV=find(elem=='VGRD');
    bandVpress=ismember(allV,pressLevs);
    bandVpress=allV(bandVpress==1);
    bandVsurf=ismember(allV,m10levs);
    bandVsurf=allV(bandVsurf==1);
    
    allP=find(elem=='PRES');
    bandPsurf=ismember(allP,m0levs);
    bandPsurf=allP(bandPsurf==1);

    tIn=readgeoraster(hourFile,Bands=bandTpress);
    rhIn=readgeoraster(hourFile,Bands=bandRHpress);
    zIn=readgeoraster(hourFile,Bands=bandZpress);
    uIn=readgeoraster(hourFile,Bands=bandUpress);
    vIn=readgeoraster(hourFile,Bands=bandVpress);
        
    %zIn=zIn./9.806;

    tps=descr(bandTpress);
    rhps=descr(bandRHpress);
    zps=descr(bandZpress);
    ups=descr(bandUpress);
    vps=descr(bandVpress);

    tp=nan(size(tps));
    zp=nan(size(zps));
    rhp=nan(size(rhps));
    up=nan(size(ups));
    vp=nan(size(vps));
    for ll=1:length(tp)
        str=char(tps(ll));
        pa=strfind(str,'[Pa]');
        tp(ll)=str2num(str(1:pa-1));

        str=char(rhps(ll));
        pa=strfind(str,'[Pa]');
        rhp(ll)=str2num(str(1:pa-1));

        str=char(zps(ll));
        pa=strfind(str,'[Pa]');
        zp(ll)=str2num(str(1:pa-1));

        str=char(ups(ll));
        pa=strfind(str,'[Pa]');
        up(ll)=str2num(str(1:pa-1));

        str=char(vps(ll));
        pa=strfind(str,'[Pa]');
        vp(ll)=str2num(str(1:pa-1));
    end

    [pCol,it]=sort(tp);
    t=tIn(:,:,it);

    [~,irh]=sort(rhp);
    rh=rhIn(:,:,irh);

    [~,iz]=sort(zp);
    z=zIn(:,:,iz);

    [~,iu]=sort(up);
    u=uIn(:,:,iu);

    [~,iv]=sort(vp);
    v=vIn(:,:,iv);

    pCol=pCol./100;
    p=repmat(pCol,1,size(t,1),size(t,2));
    p=permute(p,[2,3,1]);
    
    pHours=cat(4,pHours,p);
    tHours=cat(4,tHours,t);
    rhHours=cat(4,rhHours,rh);
    zHours=cat(4,zHours,z);
    uHours=cat(4,uHours,u);
    vHours=cat(4,vHours,v);
    
    %% Surface data
    tS=readgeoraster(hourFile,Bands=bandTsurf);
    rhS=readgeoraster(hourFile,Bands=bandRHsurf);
    uS=readgeoraster(hourFile,Bands=bandUsurf);
    vS=readgeoraster(hourFile,Bands=bandVsurf);
    pS=readgeoraster(hourFile,Bands=bandPsurf);
    
    pSurf=cat(3,pSurf,pS);
    tSurf=cat(3,tSurf,tS);
    rhSurf=cat(3,rhSurf,rhS);
    uSurf=cat(3,uSurf,uS);
    vSurf=cat(3,vSurf,vS);
    
    %% SST
    if SSTyes
        warning('No SST data available.')
    end    
end

hrrrData.Temperature=tHours;
hrrrData.rh=rhHours;
hrrrData.z=zHours;
hrrrData.u=uHours;
hrrrData.v=vHours;
hrrrData.p=pHours;
hrrrData.pSurf=pSurf;
hrrrData.tSurf=tSurf;
hrrrData.rhSurf=rhSurf;
hrrrData.uSurf=uSurf;
hrrrData.vSurf=vSurf;
hrrrData.time=inhours;
end
