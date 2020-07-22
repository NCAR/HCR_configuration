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
function [alphaTot,gammaTot,wdspd,wddir]= f_atten_layers_sfcWnd_ecmwf(soundFile,f,inlon,inlat,inalt)

%Initialize output
alphaTot=NaN;
gammaTot=NaN;
gammaTot0=NaN;
gammaTotW=NaN;

%read in lat and lon data
try
lat=ncread(soundFile,'g0_lat_0');
end
if ~exist('lat')
    try
        lat=ncread(soundFile,'g0_lat_1');
    end
end
if ~exist('lat')
    try
        lat=ncread(soundFile,'g0_lat_2');
    end
end

try
lon=ncread(soundFile,'g0_lon_1');
end
if ~exist('lon')
    try
        lon=ncread(soundFile,'g0_lon_1');
    end
end
if ~exist('lon')
    try
        lon=ncread(soundFile,'g0_lon_2');
    end
end

latInd=find(abs(lat-round(inlat,1))<0.005);
lonInd=find(abs(lon-round(inlon,1))<0.005);

temperature=squeeze(ncread(soundFile,'T_GDS0_ISBL',[lonInd,latInd,1],[1,1,inf]))-273.15;
relHum=squeeze(ncread(soundFile,'R_GDS0_ISBL',[lonInd,latInd,1],[1,1,inf]));
geopot=squeeze(ncread(soundFile,'GH_GDS0_ISBL',[lonInd,latInd,1],[1,1,inf]));

try
pressure=squeeze(ncread(soundFile,'lv_ISBL2'));
end
if ~exist('pressure')
    try
        pressure=ncread(soundFile,'lv_ISBL0');
    end
end

%surface wind speed and direction
try
sfcU=squeeze(ncread(soundFile,'10U_GDS0_SFC',[lonInd,latInd],[1,1]));
catch
    sfcU=nan;
end
try
sfcV=squeeze(ncread(soundFile,'10V_GDS0_SFC',[lonInd,latInd],[1,1]));
catch
    sfcV=nan;
end
wdspd=sqrt(sfcU^2+sfcV^2);
wddir=wrapTo360(atan2(-sfcU,-sfcV)*180/pi);

%Cat pressure, temperature, relative humidity, and altitude
ptrhalt=cat(2,pressure,temperature,relHum,geopot);

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

