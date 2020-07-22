%Reads in sounding file and calculates attenuation for each layer and sums up
%Calculates zenith one way attenuation according to the Liebe and the ITU-R method:

%Radiocommunication Sector of International Telecommunication Union. 
%Recommendation ITU-R P.676-10: Attenuation by atmospheric gases 2013.

%Liebe, H. J. (1985), An updated model for millimeter wave propagation in moist air, 
%Radio Sci., 20, 1069–1089

%Author: Ulrike Romatschke romatsch@ucar.edu
%Last modified: 20170106

%Input:
%soundFile is the path to the file containing the sounding data
%f is the radar frequency 
%varargin is an optional cut off altitude (in m above sea level). If given, the results are the
%attenuations below that altitude. The code terminates if the highest
%sounding level is below the cut off level.

%Output:
%alphaTot is the attenuation according to Liebe
%ammaTot is the attenuation according to ITU
function [alphaTot,gammaTot,wdspd,u,v]= f_atten_layers_sfcWnd_dropsonde(indir,f,intime,inalt)

%Initialize output
alphaTot=NaN;
gammaTot=NaN;
gammaTot0=NaN;
gammaTotW=NaN;
wdspd=nan;
u=nan;
v=nan;

% Find sounding file
midTime=intime(round(length(intime)));

soundFilesAll=dir([indir,'D*']);
timeDiff=[];
for ii=1:size(soundFilesAll,1)
    soundName=soundFilesAll(ii).name;
    soundTime=datetime(str2num(soundName(2:5)),str2num(soundName(6:7)),str2num(soundName(8:9)),...
        str2num(soundName(11:12)),str2num(soundName(13:14)),str2num(soundName(15:16)));
    timeDiff=cat(1,timeDiff,etime(datevec(soundTime),datevec(midTime)));
end

soundInd=find(abs(timeDiff)==min(abs(timeDiff)));
if timeDiff(soundInd)/60>10
    disp('Sounding is more than 10 minutes off!');
end

soundFile=[soundFilesAll(soundInd).folder,'/',soundFilesAll(soundInd).name];

%read in sounding data
soundDataIn=importdata(soundFile,' ',14);
soundData=soundDataIn.data;

clear soundDataIn;

%Read pressure, temperature, relative humidity, and altitude
ptrhalt=cat(2,soundData(:,5),soundData(:,6),soundData(:,8),soundData(:,14),soundData(:,11),soundData(:,9),soundData(:,10));

%Remove whole rows if there is a missing value in any of the sounding
%variables
ptrhalt(any(ptrhalt(:,1:7)==-999,2),:) = [];
ptrhalt(any(ptrhalt(:,1:7)==999,2),:) = [];
ptrhalt(any(ptrhalt(:,1:7)==99,2),:) = [];
ptrhalt(any(ptrhalt(:,1:7)==9999,2),:) = [];

if (isempty(ptrhalt)==1);
    disp('No valid sounding data.');
    return;
end

%Check if sounding is ground based or dropsonde and switch if necessary
if ptrhalt(1,4)<ptrhalt(end,4)
    ptrhalt=flipud(ptrhalt);
end

goodLevels=find(ptrhalt(:,4)<nanmean(inalt));
%Check if variables reach up to the cut off level otherwise break
if (goodLevels(1)==1);
    disp(['Plane altitude is ',num2str(nanmean(inalt)),' m but sounding ends at ',...
        num2str(ptrhalt(1,4)),' m. Results might be wrong.']);
end

%Check if cut off level equals one of the sounding levels...
if (ptrhalt(goodLevels(1),4)==nanmean(inalt));
    ptrhalt=ptrhalt(goodLevels,:);
else %...otherwise interpolate between the sounding level above and below
    %the cut off level
    pmcutoff=NaN(2,7);
    pmcutoff(1,:)=ptrhalt(goodLevels(1)-1,:);
    pmcutoff(2,:)=ptrhalt(goodLevels(1),:);
    cutoffDiff=pmcutoff(1,:)-pmcutoff(2,:);
    altperc=(nanmean(inalt)-pmcutoff(2,4))/cutoffDiff(4);
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
layer_alpha=spAlpha.*layer_depth/1000;
%sum over all layers
alphaTot=sum(layer_alpha);

%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
%multiply by layer depth
layer_gamma=spGamma.*layer_depth/1000;
%sum over all layers
gammaTot=sum(layer_gamma);

% %calculate sum of dry and wet attenuation.
% %Dry attenuation
% layer_gamma0=spGamma0.*layer_depth/1000;
% gammaTot0=sum(layer_gamma0);
% %Moist attenuation
% layer_gammaW=spGammaW.*layer_depth/1000;
% gammaTotW=sum(layer_gammaW);

%surface wind speed and direction
[sfcWndRow sfcWndCol]=find(ptrhalt(:,5)>-900 & ptrhalt(:,6)>-900 & ptrhalt(:,7)>-900);
if ptrhalt(max(sfcWndRow),4)<20;
    wdspd=ptrhalt(max(sfcWndRow),5);
    u=ptrhalt(max(sfcWndRow),6);
    v=ptrhalt(max(sfcWndRow),7);
else
    wdspd=nan;
    u=nan;
    v=nan;
end

end

