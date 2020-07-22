%Reads in sounding file and calculates attenuation for each layer and sums up
%Calculates zenith one way attenuation according to the Liebe and the ITU-R method:

%Radiocommunication Sector of International Telecommunication Union. 
%Recommendation ITU-R P.676-10: Attenuation by atmospheric gases 2013.

%Liebe, H. J. (1985), An updated model for millimeter wave propagation in moist air, 
%Radio Sci., 20, 1069–1089

%Author: Ulrike Romatschke romatsch@ucar.edu
%Last modified: 20180306

%Input:
%soundFile is the path to the file containing the sounding data
%f is the radar frequency 
%varargin is an optional cut off altitude (in m above sea level). If given, the results are the
%attenuations below that altitude. The code terminates if the highest
%sounding level is below the cut off level.

%Output:
%alphaTot is the attenuation according to Liebe
%ammaTot is the attenuation according to ITU
function [alphaTot,gammaTot,wdspd,wddir,distance]= f_atten_layers_sfcWnd_ecmwf_reanalysis(indir,f,inlon,inlat,inalt,intime)

inlon=wrapTo360(inlon);
g0=9.806;
refTime=datetime(1900,1,1,0,0,0);

%Initialize output
alphaTot=NaN;
gammaTot=NaN;
gammaTot0=NaN;
gammaTotW=NaN;
wdspd=nan;
wddir=nan;
distance=nan;

%% Pressure level data

% Find ecmwf files
roundTime=dateshift(intime, 'start', 'hour', 'nearest');
dayStr=datestr(roundTime,'yyyymmdd');

rFiles=dir([indir,'*.R.*',dayStr,'00_',dayStr,'23.nc']);
tFiles=dir([indir,'*.T.*',dayStr,'00_',dayStr,'23.nc']);
zFiles=dir([indir,'*.Z.*',dayStr,'00_',dayStr,'23.nc']);

if size(zFiles,1)==0 | size(zFiles,1)==0 | size(zFiles,1)==0
    disp('No model data found.');
    return
end

info=ncinfo([zFiles(1).folder,'/',zFiles(1).name]);

%read in lat and lon data
lonRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'longitude');
latRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'latitude');

lonInd=find(abs(lonRean-inlon)==min(abs(lonRean-inlon)));
latInd=find(abs(latRean-inlat)==min(abs(latRean-inlat)));

lonRean1=lonRean(lonInd);
latRean1=latRean(latInd);

distance=lldistkm([inlat inlon],[latRean1 lonRean1]);

ptrhalt=nan(size(tFiles,1)+1,4);

for ii=1:size(tFiles,1)
    timeRean=ncread([rFiles(ii).folder,'/',rFiles(ii).name],'time');
    timeActual=refTime+hours(timeRean);
    timeInd=find(timeActual==roundTime);
    ptrhalt(ii,1)=ncread([rFiles(ii).folder,'/',rFiles(ii).name],'level');
    ptrhalt(ii,2)=squeeze(ncread([tFiles(ii).folder,'/',tFiles(ii).name],'T',[lonInd,latInd,1,timeInd],[1,1,1,1]))-273.15;
    ptrhalt(ii,3)=squeeze(ncread([rFiles(ii).folder,'/',rFiles(ii).name],'R',[lonInd,latInd,1,timeInd],[1,1,1,1]));
    ptrhalt(ii,4)=squeeze(ncread([zFiles(ii).folder,'/',zFiles(ii).name],'Z',[lonInd,latInd,1,timeInd],[1,1,1,1]));
end

ptrhalt(:,4)=ptrhalt(:,4)./g0;

%% Surface data
% Find ecmwf files
monthStr=datestr(roundTime,'yyyymm');

dFiles=dir([indir,'*.VAR_2D.*',monthStr,'0100_*.nc']);
tsFiles=dir([indir,'*.VAR_2T.*',monthStr,'0100_*.nc']);
uFiles=dir([indir,'*.VAR_10U.*',monthStr,'0100_*.nc']);
vFiles=dir([indir,'*.VAR_10V.*',monthStr,'0100_*.nc']);
pFiles=dir([indir,'*.SP.*',monthStr,'0100_*.nc']);

if size(dFiles,1)==0 | size(tsFiles,1)==0 | size(pFiles,1)==0
    ptrhalt(end,:)=[];
else    
    info=ncinfo([dFiles.folder,'/',dFiles.name]);
    
    timeReanS=ncread([dFiles.folder,'/',dFiles.name],'time');
    timeActualS=refTime+hours(timeReanS);
    timeIndS=find(timeActualS==roundTime);
    surfPress=ncread([pFiles.folder,'/',pFiles.name],'SP',[lonInd,latInd,timeIndS],[1,1,1])./100;
    if surfPress<ptrhalt(end-1,1) & surfPress>ptrhalt(end-2,1)
        ptrhalt(end,:)=[];
        ptrhalt(end,:)=nan;
    elseif surfPress<ptrhalt(end-2,1)
        disp('Surface pressure too low.');
        ptrhalt(end,:)=[];
        return
    end
    ptrhalt(end,1)=surfPress;
    ptrhalt(end,2)=ncread([tsFiles.folder,'/',tsFiles.name],'VAR_2T',[lonInd,latInd,timeIndS],[1,1,1])-273.15;
    
    td=ncread([dFiles.folder,'/',dFiles.name],'VAR_2D',[lonInd,latInd,timeIndS],[1,1,1])-273.15;
    ptrhalt(end,3)=100*(exp((17.625*td)/(243.04+td))/exp((17.625*ptrhalt(end,2))/(243.04+ptrhalt(end,2))));
    ptrhalt(end,4)=0;
end

if size(uFiles,1)>0 & size(vFiles,1)>0
    sfcU=ncread([uFiles.folder,'/',uFiles.name],'VAR_10U',[lonInd,latInd,timeIndS],[1,1,1]);
    sfcV=ncread([vFiles.folder,'/',vFiles.name],'VAR_10V',[lonInd,latInd,timeIndS],[1,1,1]);
    
    wdspd=sqrt(sfcU^2+sfcV^2);
    wddir=wrapTo360(atan2(-sfcU,-sfcV)*180/pi);
end

%% Attenuation

%Check if sounding is ground based or dropsonde and switch if necessary
if ptrhalt(1,4)<ptrhalt(end,4)
    ptrhalt=flipud(ptrhalt);
end

%Cut off above flight level
cutoff=inalt;
goodLevels=find(ptrhalt(:,4)<cutoff);
%Check if variables reach up to the cut off level otherwise break
if (goodLevels(1)==1);
    disp('Sounding data end below the cut off level.');
    return;
end

%Check if cut off level equals one of the sounding levels...
if (ptrhalt(goodLevels(1),4)==cutoff);
    ptrhalt=ptrhalt(goodLevels,:);
else %...otherwise interpolate between the sounding level above and below
    %the cut off level
    pmcutoff=NaN(2,4);
    pmcutoff(1,:)=ptrhalt(goodLevels(1)-1,:);
    pmcutoff(2,:)=ptrhalt(goodLevels(1),:);
    cutoffDiff=pmcutoff(1,:)-pmcutoff(2,:);
    altperc=(cutoff-pmcutoff(2,4))/cutoffDiff(4);
    cutoffVals=pmcutoff(2,:)+cutoffDiff*altperc;
    ptrhalt=cat(1,cutoffVals,ptrhalt(goodLevels,:));
end

%calculate mean of the variables for each layer between the sounding levels
layer_p=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,1),ptrhalt(2:size(ptrhalt,1),1)),2);
layer_t=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,2),ptrhalt(2:size(ptrhalt,1),2)),2);
layer_rh=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,3),ptrhalt(2:size(ptrhalt,1),3)),2);
layer_depth=abs(ptrhalt(1:size(ptrhalt,1)-1,4)-ptrhalt(2:size(ptrhalt,1),4));
%saturation pressure after WMO
%layer_es=6.112.*exp((17.62.*layer_t)./(layer_t+243.12));
%saturation pressure after Liebe 1981
layer_es=2.409.*(300./(layer_t+273.15)).^5.*10.^(10-9.834.*(300./(layer_t+273.15))).*10;

%run Liebe code which gives attenuation per kilometer
[spAlpha]=f_atten_Liebe(f,layer_p,layer_t,layer_es,layer_rh);
%multiply by layer depth
layer_alpha=spAlpha.*double(layer_depth)/1000;
%sum over all layers
alphaTot=sum(layer_alpha);

%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
%multiply by layer depth
layer_gamma=spGamma.*double(layer_depth)/1000;
%sum over all layers
gammaTot=sum(layer_gamma);

% %calculate sum of dry and wet attenuation.
% %Dry attenuation
% layer_gamma0=spGamma0.*layer_depth/1000;
% gammaTot0=sum(layer_gamma0);
% %Moist attenuation
% layer_gammaW=spGammaW.*layer_depth/1000;
% gammaTotW=sum(layer_gammaW);


end

