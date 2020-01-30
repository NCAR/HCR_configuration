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
function [alphaTot,gammaTot,gammaTot0,gammaTotW]= f_atten_layers(soundFile,f,varargin)

%Initialize output
alphaTot=NaN;
gammaTot=NaN;
gammaTot0=NaN;
gammaTotW=NaN;

%read in sounding data
soundDataIn=importdata(soundFile,' ',14);
soundData=soundDataIn.data;

clear soundDataIn;

%Read pressure, temperature, relative humidity, and altitude
ptrhalt=cat(2,soundData(:,5),soundData(:,6),soundData(:,8),soundData(:,14));

%Remove whole rows if there is a missing value in any of the sounding
%variables
ptrhalt(any(ptrhalt==-999,2),:) = [];

if (isempty(ptrhalt)==1);
    disp('No valid sounding data.');
    return;
end

%Check if sounding is ground based or dropsonde and switch if necessary
if ptrhalt(1,4)<ptrhalt(end,4)
    ptrhalt=flipud(ptrhalt);
end

%Check if cut off value is given
if (any(size(varargin))==1);
    cutoff=varargin{1};
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

%calculate sum of dry and wet attenuation.
%Dry attenuation
layer_gamma0=spGamma0.*layer_depth/1000;
gammaTot0=sum(layer_gamma0);
%Moist attenuation
layer_gammaW=spGammaW.*layer_depth/1000;
gammaTotW=sum(layer_gammaW);

