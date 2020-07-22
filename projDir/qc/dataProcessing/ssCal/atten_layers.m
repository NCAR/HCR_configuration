%Reads in sounding file and calculates attenuation for each layer and sums up
%Calculates attenuation according to the Liebe and the ITU-R method:

%Radiocommunication Sector of International Telecommunication Union. 
%Recommendation ITU-R P.676-10: Attenuation by atmospheric gases 2013.

%Liebe, H. J. (1985), An updated model for millimeter wave propagation in moist air, 
%Radio Sci., 20, 1069–1089

%Author: Ulrike Romatschke romatsch@ucar.edu
%Last modified: 20170106

clear all;
close all;

%Directory with functions (if not in same directory)
addpath('/h/eol/romatsch/matlab/radar/functions/');

%Sounding data file
soundFile='/h/eol/romatsch/data/radar/soundings/D20150807_212000_P.QC.eol';

f = 61; %94.4; %Radar frequency in GHz
angle=0; %beamangle from vertical in degree

%Import sounding data
soundDataIn=importdata(soundFile,' ',13);
soundData=soundDataIn.data;

clear soundDataIn;

ptrhalt=cat(2,soundData(:,5),soundData(:,6),soundData(:,8),soundData(:,14));

ptrhalt(any(ptrhalt==-999,2),:) = [];

layer_p=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,1),ptrhalt(2:size(ptrhalt,1),1)),2);
layer_t=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,2),ptrhalt(2:size(ptrhalt,1),2)),2);
layer_rh=mean(cat(2,ptrhalt(1:size(ptrhalt,1)-1,3),ptrhalt(2:size(ptrhalt,1),3)),2);
layer_depth=abs(ptrhalt(1:size(ptrhalt,1)-1,4)-ptrhalt(2:size(ptrhalt,1),4));
%saturation pressure after WMO
%layer_es=6.112.*exp((17.62.*layer_t)./(layer_t+243.12));
%saturation pressure after Liebe 1981
layer_es=2.409.*(300./(layer_t+273.15)).^5.*10.^(10-9.834.*(300./(layer_t+273.15))).*10;

%Liebe code
[spAlpha]=f_atten_Liebe(f,layer_p,layer_t,layer_es,layer_rh);
layer_alpha=spAlpha.*layer_depth/1000;
alphaTot=sum(layer_alpha);
alphaAngle=alphaTot/cos(degtorad(angle));

%ITU-R code
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
%Total attenuation
layer_gamma=spGamma.*layer_depth/1000;
gammaTot=sum(layer_gamma);
gammaAngle=gammaTot/cos(degtorad(angle));
%Dry attenuation
layer_gamma0=spGamma0.*layer_depth/1000;
gammaTot0=sum(layer_gamma0);
gammaAngle0=gammaTot0/cos(degtorad(angle));
%Moist attenuation
layer_gammaW=spGammaW.*layer_depth/1000;
gammaTotW=sum(layer_gammaW);
gammaAngleW=gammaTotW/cos(degtorad(angle));

outputLiebe=['Attenuation with Liebe method is ',num2str(round(alphaAngle,2)),' dB.'];
outputITUR=['Attenuation with ITU-R method is ',num2str(round(gammaAngle,2)),' dB.'];

disp(outputLiebe);
disp(outputITUR);