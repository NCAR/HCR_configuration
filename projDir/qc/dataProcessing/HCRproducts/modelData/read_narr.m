function narrData = read_narr(indir,startTime,endTime,SSTyes)
% Find the right time span for era5 data and read them in

refTime=datetime(1900,1,1,0,0,0);

%Initialize output
narrData=[];

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


% Get intime hours
firstHour=hour(startTime);
threeHour=floor(firstHour/3)*3;

startHour=datetime(year(startTime),month(startTime),day(startTime),threeHour,0,0);
endHour=datetime(year(endTime),month(endTime),day(endTime),hour(endTime)+3,0,0);

inhours=startHour:hours(3):endHour;

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

modelFiles3D=dir([indir,'NARR3D*.nc']);
fileNames3D=cell2mat({modelFiles3D.name}');

startFileTimes3D=datetime(str2num(fileNames3D(:,8:11)),str2num(fileNames3D(:,12:13)),str2num(fileNames3D(:,15:16)));

startFileInds3D=min(find(startFileTimes3D>startHour))-1;
endFileInds3D=max(find(startFileTimes3D<=endHour));
if isempty(startFileInds3D) & ~isempty(endFileInds3D)
    startFileInds3D=endFileInds3D;
end

rightFiles3D=fileNames3D(startFileInds3D:endFileInds3D,:);

if size(rightFiles3D,1)==0
    disp('No pressure level model data found.');
    return
end

modelFilesSurf=dir([indir,'NARRflx*.nc']);
fileNamesSurf=cell2mat({modelFilesSurf.name}');

startFileTimesSurf=datetime(str2num(fileNamesSurf(:,9:12)),str2num(fileNamesSurf(:,13:14)),str2num(fileNamesSurf(:,16:17)));

startFileIndsSurf=min(find(startFileTimesSurf>startHour))-1;
endFileIndsSurf=max(find(startFileTimesSurf<=endHour));
if isempty(startFileIndsSurf) & ~isempty(endFileIndsSurf)
    startFileIndsSurf=endFileIndsSurf;
end

rightFilesSurf=fileNamesSurf(startFileIndsSurf:endFileIndsSurf,:);

if size(rightFilesSurf,1)==0
    disp('No surface model data found.');
    return
end

% Read lon and lat

lonRean=ncread([indir,rightFiles3D(1,:)],'gridlon_221');
lonRean=wrapTo360(lonRean);
latRean=ncread([indir,rightFiles3D(1,:)],'gridlat_221');
levRean=ncread([indir,rightFiles3D(1,:)],'lv_ISBL1');

%%%%%%%%%%%%%%%%%%% start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:size(rightFiles3D,1)
    disp(rightFiles3D(ii,:));
    timeReanIn=ncread([indir,rightFiles3D(ii,:)],'initial_time0_hours');
    timeRean=datetime(1800,1,1)+hours(timeReanIn);
    
    [~,timeIndsNeeded,~]=intersect(timeRean,inhours);
    
    tIn=ncread([indir,rightFiles3D(ii,:)],'TMP_221_ISBL');
    t=tIn(:,:,:,timeIndsNeeded);
    pIn=repmat(levRean,1,size(t,2),size(t,1),size(t,4));
    p=double(permute(pIn,[3,2,1,4]));
    shIn=ncread([indir,rightFiles3D(ii,:)],'SPF_H_221_ISBL');
    sh=shIn(:,:,:,timeIndsNeeded);
    rh=26.3.*p.*sh./(exp((17.67.*(t-273.15))./(t-29.65)));
    t=t-273.15;
    zIn=ncread([indir,rightFiles3D(ii,:)],'HGT_221_ISBL');
    z=zIn(:,:,:,timeIndsNeeded);
    
    pHours=cat(4,pHours,p);
    tHours=cat(4,tHours,t);
    rhHours=cat(4,rhHours,rh);
    zHours=cat(4,zHours,z);
end
    
levReanS=ncread([indir,rightFilesSurf(1,:)],'lv_HTGL1');

for ii=1:size(rightFilesSurf,1)
    disp(rightFilesSurf(ii,:));
    timeReanIn=ncread([indir,rightFilesSurf(ii,:)],'initial_time0_hours');
    timeRean=datetime(1800,1,1)+hours(timeReanIn);
    
    [~,timeIndsNeeded,~]=intersect(timeRean,inhours);
    
    tsIn=ncread([indir,rightFilesSurf(ii,:)],'TMP_221_HTGL');
    tS=squeeze(tsIn(:,:,1,timeIndsNeeded));
    tS=tS-273.15;
    psIn=ncread([indir,rightFilesSurf(ii,:)],'PRES_221_HTGL');
    pS=squeeze(psIn(:,:,1,timeIndsNeeded)./100);
    tdsIn=ncread([indir,rightFilesSurf(ii,:)],'DPT_221_HTGL');
    tdS=squeeze(tdsIn(:,:,timeIndsNeeded));
    tdS=tdS-273.15;
    rhS=100*(exp((17.625*tdS)./(243.04+tdS))./exp((17.625*tS(:,:,end))./(243.04+tS(:,:,end))));
    usIn=ncread([indir,rightFilesSurf(ii,:)],'U_GRD_221_HTGL');
    sfcU=squeeze(usIn(:,:,timeIndsNeeded));
    vsIn=ncread([indir,rightFilesSurf(ii,:)],'V_GRD_221_HTGL');
    sfcV=squeeze(vsIn(:,:,timeIndsNeeded));
        
    pSurf=cat(3,pSurf,pS);
    tSurf=cat(3,tSurf,tS);
    rhSurf=cat(3,rhSurf,rhS);
    uSurf=cat(3,uSurf,sfcU);
    vSurf=cat(3,vSurf,sfcV);
    
    %% SST
    if SSTyes
        disp('Getting SST is not set up yet.');
    end
end

narrData.lat=latRean;
narrData.lon=lonRean;
narrData.Temperature=tHours;
narrData.rh=rhHours;
narrData.z=zHours;
narrData.p=pHours;
narrData.pSurf=pSurf;
narrData.tSurf=tSurf;
narrData.rhSurf=rhSurf;
narrData.uSurf=uSurf;
narrData.vSurf=vSurf;
narrData.time=inhours;
if SSTyes
    narrData.sstSurf=sstSurf;
end

end

