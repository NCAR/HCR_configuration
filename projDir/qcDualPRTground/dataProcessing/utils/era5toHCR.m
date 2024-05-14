%Interpolate era5 data to HCR track

function [era5fields]= era5toHCR(indir,inlon,inlat,intime,SSTyes)
if size(inlon,1)~=1
    inlon=inlon';
    inlat=inlat';
    intime=intime';
end

inlon=wrapTo360(inlon);
refTime=datetime(1900,1,1,0,0,0);

%Initialize output
era5fields=[];

% Get intime hours
inhoursAll=datetime(year(intime),month(intime),day(intime),hour(intime),0,0);
inhoursTemp=unique(inhoursAll);
inhoursLast=inhoursTemp(end)+hours(1);

ll=1;
inhours=inhoursTemp(1);
while inhours(end)~=inhoursLast
    inhours=[inhours inhoursTemp(1)+hours(ll)];
    ll=ll+1;
end

pHours={};
tHours={};
rhHours={};
zHours={};
uHours={};
vHours={};
if SSTyes
    sstHours={};
end

for kk=1:length(inhours) %Loop through all hours
    %% Pressure level data
    % Find ecmwf files
    roundTime=inhours(kk);
    dayStr=datestr(roundTime,'yyyymmdd');
    
    rFiles=dir([indir,'*.R.*',dayStr,'00_',dayStr,'23.nc']);
    tFiles=dir([indir,'*.T.*',dayStr,'00_',dayStr,'23.nc']);
    zFiles=dir([indir,'*.Z.*',dayStr,'00_',dayStr,'23.nc']);
    
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
    
    t=nan(length(lonRean),length(latRean),length(rFiles)+1);
    rh=nan(length(lonRean),length(latRean),length(rFiles)+1);
    z=nan(length(lonRean),length(latRean),length(rFiles)+1);
    p=nan(length(lonRean),length(latRean),length(rFiles)+1);
    
    for jj=1:size(tFiles,1)
        p(:,:,jj)=fliplr(ncread([rFiles(jj).folder,'/',rFiles(jj).name],'level'));
        t(:,:,jj)=fliplr(squeeze(ncread([tFiles(jj).folder,'/',tFiles(jj).name],'T',[1,1,1,timeInd],[inf,inf,inf,1]))-273.15);
        rh(:,:,jj)=fliplr(squeeze(ncread([rFiles(jj).folder,'/',rFiles(jj).name],'R',[1,1,1,timeInd],[inf,inf,inf,1])));
        z(:,:,jj)=fliplr(squeeze(ncread([zFiles(jj).folder,'/',zFiles(jj).name],'Z',[1,1,1,timeInd],[inf,inf,inf,1])));
    end
    
    z=z./9.806;
    %% Surface data
    % Find ecmwf files
    monthStr=datestr(roundTime,'yyyymm');
    
    dFiles=dir([indir,'*.VAR_2D.*',monthStr,'0100_*.nc']);
    tsFiles=dir([indir,'*.VAR_2T.*',monthStr,'0100_*.nc']);
    uFiles=dir([indir,'*.VAR_10U.*',monthStr,'0100_*.nc']);
    vFiles=dir([indir,'*.VAR_10V.*',monthStr,'0100_*.nc']);
    pFiles=dir([indir,'*.SP.*',monthStr,'0100_*.nc']);
    
    info=ncinfo([dFiles.folder,'/',dFiles.name]);
    
    timeReanS=ncread([dFiles.folder,'/',dFiles.name],'time');
    timeActualS=refTime+hours(timeReanS);
    timeIndS=find(timeActualS==roundTime);
    p(:,:,end)=fliplr(ncread([pFiles.folder,'/',pFiles.name],'SP',[1,1,timeIndS],[inf,inf,1])./100);
    t(:,:,end)=fliplr(ncread([tsFiles.folder,'/',tsFiles.name],'VAR_2T',[1,1,timeIndS],[inf,inf,1])-273.15);
    td=fliplr(ncread([dFiles.folder,'/',dFiles.name],'VAR_2D',[1,1,timeIndS],[inf,inf,1])-273.15);
    rh(:,:,end)=100*(exp((17.625*td)./(243.04+td))./exp((17.625*t(:,:,end))./(243.04+t(:,:,end))));
    z(:,:,end)=0;
    
    % Sort out data that is below sea surface
    pSurf=p(:,:,end);
    tSurf=t(:,:,end);
    rhSurf=rh(:,:,end);
    zSurf=z(:,:,end);
    
    % If only lowest level above surface level
    for ii=1:length(rFiles)
        %lowinds=find(p(:,:,end)<p(:,:,end-1));
        pTemp=p(:,:,ii);
        lowinds=find(pSurf<pTemp);
        if ii>1
            naninds=intersect(lowinds,lowindsOld);
        else
            naninds=[];
        end
        pTemp(lowinds)=pSurf(lowinds);
        pTemp(naninds)=nan;
        p(:,:,ii)=pTemp;
        
        tTemp=t(:,:,ii);
        tTemp(lowinds)=tSurf(lowinds);
        tTemp(naninds)=nan;
        t(:,:,ii)=tTemp;
        
        rhTemp=rh(:,:,ii);
        rhTemp(lowinds)=rhSurf(lowinds);
        rhTemp(naninds)=nan;
        rh(:,:,ii)=rhTemp;
        
        zTemp=z(:,:,ii);
        zTemp(lowinds)=zSurf(lowinds);
        zTemp(naninds)=nan;
        zNo0=zTemp;
        zeroInds=find(zTemp==0);
        zNo0(zeroInds)=nan;
        zInterp = fillmissing(zNo0,'linear','EndValues',0);
        zTemp(zeroInds)=zInterp(zeroInds);
        
        z(:,:,ii)=zTemp;
        
        if ii>1
            indsTemp=cat(1,lowinds,lowindsOld);
            lowinds=unique(indsTemp);
        end
        
        lowindsOld=lowinds;
    end
    
    if size(uFiles,1)>0 & size(vFiles,1)>0
        sfcU=fliplr(ncread([uFiles.folder,'/',uFiles.name],'VAR_10U',[1,1,timeIndS],[inf,inf,1]));
        sfcV=fliplr(ncread([vFiles.folder,'/',vFiles.name],'VAR_10V',[1,1,timeIndS],[inf,inf,1]));
    end
    
    %% SST
    if SSTyes
        sstFiles=dir([indir,'*.SSTK.*',monthStr,'0100_*.nc']);
        sst=fliplr(ncread([sstFiles.folder,'/',sstFiles.name],'SSTK',[1,1,timeIndS],[inf,inf,1])-273.15);
        sstHours{end+1}=sst;
    end        
        
    pHours{end+1}=p;
    tHours{end+1}=t;
    rhHours{end+1}=rh;
    zHours{end+1}=z;
    uHours{end+1}=sfcU;
    vHours{end+1}=sfcV;
end

lon2d=repmat(lonRean,1,length(latRean));
lat2d=repmat(fliplr(latRean'),length(lonRean),1);

pHCR=nan(length(inlon),size(p,3));
tHCR=nan(length(inlon),size(p,3));
rhHCR=nan(length(inlon),size(p,3));
zHCR=nan(length(inlon),size(p,3));

%% Interpolate era5 to sub intervals
% Create time vector
MinInt=1;

ll=1;
timesub=inhours(1);
while timesub(end)~=inhours(end)
    timesub=[timesub inhours(1)+minutes(ll*MinInt)];
    ll=ll+1;
end

psub={};
tsub={};
rhsub={};
zsub={};
usub={};
vsub={};
if SSTyes
    sstsub={};
end

perHour=60/MinInt;
perHour1=perHour-1;
x=repmat([0:1:perHour1],1,length(inhours));

for kk=1:length(timesub)
    startHour=floor((kk-1)/perHour)+1;
    endHour=floor((kk-1)/perHour)+2;
    
    if x(kk)==0
        psub{end+1}=pHours{1,startHour};
        tsub{end+1}=tHours{1,startHour};
        rhsub{end+1}=rhHours{1,startHour};
        zsub{end+1}=zHours{1,startHour};
        usub{end+1}=uHours{1,startHour};
        vsub{end+1}=vHours{1,startHour};
        if SSTyes
            sstsub{end+1}=sstHours{1,startHour};
        end
    else
        psub{end+1}=(pHours{1,startHour}.*(perHour1-x(kk))+pHours{1,endHour}.*x(kk))./perHour1;
        tsub{end+1}=(tHours{1,startHour}.*(perHour1-x(kk))+tHours{1,endHour}.*x(kk))./perHour1;
        rhsub{end+1}=(rhHours{1,startHour}.*(perHour1-x(kk))+rhHours{1,endHour}.*x(kk))./perHour1;
        zsub{end+1}=(zHours{1,startHour}.*(perHour1-x(kk))+zHours{1,endHour}.*x(kk))./perHour1;
        usub{end+1}=(uHours{1,startHour}.*(perHour1-x(kk))+uHours{1,endHour}.*x(kk))./perHour1;
        vsub{end+1}=(vHours{1,startHour}.*(perHour1-x(kk))+vHours{1,endHour}.*x(kk))./perHour1;
        if SSTyes
            sstsub{end+1}=(sstHours{1,startHour}.*(perHour1-x(kk))+sstHours{1,endHour}.*x(kk))./perHour1;
        end
    end
end

%% Find closest time
hcrTimeMat=repmat(intime,length(timesub),1);
timesubMat=repmat(timesub,length(intime),1)';
timeDiff=etime(datevec(hcrTimeMat),datevec(timesubMat));
timeDiffMat=abs(reshape(timeDiff,size(timesubMat)));

[minDiff minInd]=min(timeDiffMat,[],1);
minIndU=unique(minInd);

for kk=1:length(minIndU)
    inInds=find(minInd==minIndU(kk));
    
    lonInds=inlon(inInds);
    latInds=inlat(inInds);
    
    for ii=1:size(p,3)
        Fp=griddedInterpolant(lon2d,lat2d,psub{1,minIndU(kk)}(:,:,ii));
        pHCR(inInds,ii)=Fp(lonInds,latInds);
        Ft=griddedInterpolant(lon2d,lat2d,tsub{1,minIndU(kk)}(:,:,ii));
        tHCR(inInds,ii)=Ft(lonInds,latInds);
        Frh=griddedInterpolant(lon2d,lat2d,rhsub{1,minIndU(kk)}(:,:,ii));
        rhHCR(inInds,ii)=Frh(lonInds,latInds);
        Fz=griddedInterpolant(lon2d,lat2d,zsub{1,minIndU(kk)}(:,:,ii));
        zHCR(inInds,ii)=Fz(lonInds,latInds);
    end
    
    Fu=griddedInterpolant(lon2d,lat2d,usub{1,minIndU(kk)});
    uHCR(inInds)=Fu(lonInds,latInds);
    Fv=griddedInterpolant(lon2d,lat2d,vsub{1,minIndU(kk)});
    vHCR(inInds)=Fv(lonInds,latInds);
    if SSTyes
        Fsst=griddedInterpolant(lon2d,lat2d,sstsub{1,minIndU(kk)});
        sstHCR(inInds)=Fsst(lonInds,latInds);
    end
end

era5fields.p=pHCR;
era5fields.t=tHCR;
era5fields.rh=rhHCR;
era5fields.z=zHCR;
era5fields.u=uHCR';
era5fields.v=vHCR';
if SSTyes
    era5fields.sst=sstHCR';
end
% close all
% figure
% fig1=surf(lon2d,lat2d,rhsub{1,1}(:,:,6));
% %fig1=surf(lon2d,lat2d,usub{1,1}(:,:));
% fig1.EdgeColor='none';
% view(2);
% caxis([0 110]);
% hold on
% plot(inlon,inlat,'k','linewidth',2);
% ax = gca;
% ax.SortMethod = 'childorder';
% colorbar
% 
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
% 
% set(gcf,'PaperPositionMode','auto')
% print(['/h/eol/romatsch/hcrCalib/notes/rhum_lonlat'],'-dpng','-r0');

% figure
% if size(intime,2)~=1
%     intime=intime';
% end
% fig1=surf(repmat(intime,1,size(pHCR,2)),zHCR./1000,rhHCR);
% %fig1=surf(repmat(intime,1,size(pHCR,2)),zHCR./1000,vHCR);
% fig1.EdgeColor='none';
% caxis([0 110]);
% xlim([intime(1),intime(end)]);
% ylabel('Altitude (km)');
% view(2);
% colorbar
% 
% set(gcf,'PaperPositionMode','auto')
% print(['/h/eol/romatsch/hcrCalib/notes/rhum_cross'],'-dpng','-r0');
end

