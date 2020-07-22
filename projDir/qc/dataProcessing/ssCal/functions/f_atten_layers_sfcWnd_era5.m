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
function [alphaTot,gammaTot,wdspd,u,v,returnInds]= f_atten_layers_sfcWnd_era5(indir,f,inlon,inlat,inalt,intime)

inlon=wrapTo360(inlon);
g0=9.806;
refTime=datetime(1900,1,1,0,0,0);

%Initialize output
alphaTot=nan(size(inlon));
gammaTot=nan(size(inlon));
gammaTot0=nan(size(inlon));
gammaTotW=nan(size(inlon));
wdspd=nan(size(inlon));
u=nan(size(inlon));
v=nan(size(inlon));
distance=nan(size(inlon));
ulonlatInd=nan;

%% Pressure level data

% Find ecmwf files
roundTime=dateshift(intime(end), 'start', 'hour', 'nearest');
dayStr=datestr(roundTime,'yyyymmdd');

rFiles=dir([indir,'*.R.*',dayStr,'00_',dayStr,'23.nc']);
tFiles=dir([indir,'*.T.*',dayStr,'00_',dayStr,'23.nc']);
zFiles=dir([indir,'*.Z.*',dayStr,'00_',dayStr,'23.nc']);

if size(zFiles,1)==0 | size(zFiles,1)==0 | size(zFiles,1)==0
    disp('No model data found.');
    return
end

%read in lat and lon data
lonRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'longitude');
latRean=ncread([rFiles(1).folder,'/',rFiles(1).name],'latitude');

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
returnInds={ulonlatInd,ic};

for ii=1:length(ia)
    lonRean1=lonRean(ulonlatInd(ii,1));
    latRean1=latRean(ulonlatInd(ii,2));
    
    distance(ic==ii)=lldistkm([inlat inlon],[latRean1 lonRean1]);
    
    ptrhalt=nan(size(tFiles,1)+1,4);
    
    for jj=1:size(tFiles,1)
        timeRean=ncread([rFiles(jj).folder,'/',rFiles(jj).name],'time');
        timeActual=refTime+hours(timeRean);
        timeInd=find(timeActual==roundTime);
        ptrhalt(jj,1)=ncread([rFiles(jj).folder,'/',rFiles(jj).name],'level');
        ptrhalt(jj,2)=squeeze(ncread([tFiles(jj).folder,'/',tFiles(jj).name],'T',[ulonlatInd(ii,1),ulonlatInd(ii,2),1,timeInd],[1,1,1,1]))-273.15;
        ptrhalt(jj,3)=squeeze(ncread([rFiles(jj).folder,'/',rFiles(jj).name],'R',[ulonlatInd(ii,1),ulonlatInd(ii,2),1,timeInd],[1,1,1,1]));
        ptrhalt(jj,4)=squeeze(ncread([zFiles(jj).folder,'/',zFiles(jj).name],'Z',[ulonlatInd(ii,1),ulonlatInd(ii,2),1,timeInd],[1,1,1,1]));
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
        surfPress=ncread([pFiles.folder,'/',pFiles.name],'SP',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1])./100;
        if surfPress<ptrhalt(end-1,1) & surfPress>ptrhalt(end-2,1)
            ptrhalt(end,:)=[];
            ptrhalt(end,:)=nan;
        elseif surfPress<ptrhalt(end-2,1)
            disp('Surface pressure too low.');
            ptrhalt(end,:)=[];
            return
        end
        ptrhalt(end,1)=surfPress;
        ptrhalt(end,2)=ncread([tsFiles.folder,'/',tsFiles.name],'VAR_2T',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1])-273.15;
        
        td=ncread([dFiles.folder,'/',dFiles.name],'VAR_2D',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1])-273.15;
        ptrhalt(end,3)=100*(exp((17.625*td)/(243.04+td))/exp((17.625*ptrhalt(end,2))/(243.04+ptrhalt(end,2))));
        ptrhalt(end,4)=0;
    end
    
    if size(uFiles,1)>0 & size(vFiles,1)>0
        sfcU=ncread([uFiles.folder,'/',uFiles.name],'VAR_10U',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1]);
        sfcV=ncread([vFiles.folder,'/',vFiles.name],'VAR_10V',[ulonlatInd(ii,1),ulonlatInd(ii,2),timeIndS],[1,1,1]);
        
        wdspd(ic==ii)=sqrt(sfcU^2+sfcV^2);
        u(ic==ii)=sfcU;
        v(ic==ii)=sfcV;
    end
    
    %% Attenuation
    
    %Check if sounding is ground based or dropsonde and switch if necessary
    if ptrhalt(1,4)<ptrhalt(end,4)
        ptrhalt=flipud(ptrhalt);
    end
    
    %Cut off above flight level
    cutoff=nanmean(inalt);
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
    alphaTot(ic==ii)=sum(layer_alpha);
    
    %run ITU-R code which gives total, dry, and wet, attenuation per kilometer
    [spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
    %multiply by layer depth
    layer_gamma=spGamma.*double(layer_depth)/1000;
    %sum over all layers
    gammaTot(ic==ii)=sum(layer_gamma);
    
    % %calculate sum of dry and wet attenuation.
    % %Dry attenuation
    % layer_gamma0=spGamma0.*layer_depth/1000;
    % gammaTot0=sum(layer_gamma0);
    % %Moist attenuation
    % layer_gammaW=spGammaW.*layer_depth/1000;
    % gammaTotW=sum(layer_gammaW);
end

end

