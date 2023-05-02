%
%Radiocommunication Sector of International Telecommunication Union. 
%Recommendation ITU-R P.676-10: Attenuation by atmospheric gases 2013.

%Liebe, H. J. (1985), An updated model for millimeter wave propagation in moist air, 
%Radio Sci., 20, 1069–1089

%f is the radar frequency
%varargin is an optional cut off altitude (in m above sea level). If given, the results are the
%attenuations below that altitude. The code terminates if the highest
%sounding level is below the cut off level.

%Output:
%alphaTot is the attenuation according to Liebe
%gammaTot is the attenuation according to ITU
function [ituTot,ituCloud,layer_itu,layer_ituC]= get_gas_atten(data)

f=data.frq/1e+9;

if ~isfield(data,'dbzMasked')
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
end

% Layer depth
layer_depth=abs(diff(data.asl,1,1));
layer_depth=cat(1,layer_depth,layer_depth(end,:));

% Set rh in cloud to 100 %
RHcloud=data.RH;
RHcloud(~isnan(data.dbzMasked))=100;

%saturation pressure after WMO
%layer_es=6.112.*exp((17.62.*layer_t)./(layer_t+243.12));
%saturation pressure after Liebe 1981
layer_es=2.409.*(300./(data.TEMP+273.15)).^5.*10.^(10-9.834.*(300./(data.TEMP+273.15))).*10;

%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,data.PRESS,data.TEMP,data.RH,layer_es);
%multiply by layer depth
layer_itu=spGamma.*layer_depth./1000;
%sum over all layers
ituTot=sum(layer_itu,1,'omitnan')';

% Calculate attenuation with cloud
%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGammaC, spGamma0C, spGammaWC]=f_atten_ITUR(f,data.PRESS,data.TEMP,RHcloud,layer_es);
%multiply by layer depth
layer_ituC=spGammaC.*layer_depth./1000;
%sum over all layers
ituCloud=sum(layer_ituC,1,'omitnan')';

end



