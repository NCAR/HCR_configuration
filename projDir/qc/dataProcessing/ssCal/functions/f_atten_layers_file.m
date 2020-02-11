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
function [alphaTot,gammaTot,wdspd,u,v,SST]= f_atten_layers_file(f,intime,SSTyes,modeldata)

%Initialize output
alphaTot=nan(size(intime));
gammaTot=nan(size(intime));
wdspd=nan(size(intime));
u=nan(size(intime));
v=nan(size(intime));
SST=nan(size(intime));

wdspd(:)=sqrt(modeldata.u.^2+modeldata.v.^2);
u(:)=modeldata.u;
v(:)=modeldata.v;
if SSTyes
    SST(:)=modeldata.sst;
end

modFields=fields(modeldata);
for ii=1:length(modFields)
    modeldata.(modFields{ii})=modeldata.(modFields{ii})';
end

%% Attenuation
ptrhalt=cat(3,modeldata.p,modeldata.temp,modeldata.rh,modeldata.rangeMat);

%calculate mean of the variables for each layer between the sounding levels
layer_p=(ptrhalt(1:size(ptrhalt,1)-1,:,1)+ptrhalt(2:size(ptrhalt,1),:,1))./2;
layer_t=(ptrhalt(1:size(ptrhalt,1)-1,:,2)+ptrhalt(2:size(ptrhalt,1),:,2))./2;
layer_rh=(ptrhalt(1:size(ptrhalt,1)-1,:,3)+ptrhalt(2:size(ptrhalt,1),:,3))./2;
layer_depth=-ptrhalt(1:size(ptrhalt,1)-1,:,4)+ptrhalt(2:size(ptrhalt,1),:,4);

%saturation pressure after WMO
%layer_es=6.112.*exp((17.62.*layer_t)./(layer_t+243.12));
%saturation pressure after Liebe 1981
layer_es=2.409.*(300./(layer_t+273.15)).^5.*10.^(10-9.834.*(300./(layer_t+273.15))).*10;

%run Liebe code which gives attenuation per kilometer
[spAlpha]=f_atten_Liebe(f,layer_p,layer_t,layer_es,layer_rh);
%multiply by layer depth
layer_alpha=spAlpha.*layer_depth./1000;
% Remove data where range is negative
layer_alpha(modeldata.rangeMat(1:end-1,:)<0)=nan;
%sum over all layers
alphaTot1=sum(layer_alpha,1,'omitnan')';

%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
%multiply by layer depth
layer_gamma=spGamma.*layer_depth./1000;
% Remove data where range is negative
layer_gamma(modeldata.rangeMat(1:end-1,:)<0)=nan;
%sum over all layers
gammaTot1=sum(layer_gamma,1,'omitnan')';

%In principle, this code already computes the total attenuation but because
%it will be adjusted for angle in f_determine_bias, we artificially put it on
%a vertical path. The 5 instead of a 10 is to account for two way
%attenuation in f_determine_bias.
modeldata.elev(modeldata.elev>90)=nan;
alphaTot=alphaTot1+5*log10(cosd(modeldata.elev'));
gammaTot=gammaTot1+5*log10(cosd(modeldata.elev'));
end



