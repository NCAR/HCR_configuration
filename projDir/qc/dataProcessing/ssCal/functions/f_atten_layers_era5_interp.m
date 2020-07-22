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
function [alphaTot,gammaTot,wdspd,u,v,SST]= f_atten_layers_sfcWnd_era5(indir,f,inlon,inlat,inalt,intime,SSTyes)

%Initialize output
alphaTot=nan(size(inlon));
gammaTot=nan(size(inlon));
wdspd=nan(size(inlon));
u=nan(size(inlon));
v=nan(size(inlon));
SST=nan(size(inlon));

era5data=era5toHCR(indir,inlon,inlat,intime,SSTyes);

wdspd(:)=sqrt(era5data.u.^2+era5data.v.^2);
u(:)=era5data.u;
v(:)=era5data.v;
if SSTyes
    SST(:)=era5data.sst;
end
    
%% Attenuation
ptrhalt=cat(3,era5data.p,era5data.t,era5data.rh,era5data.z);
ptrhalt=permute(ptrhalt,[2,1,3]);
%Check if sounding is ground based or dropsonde and switch if necessary
if ptrhalt(1,1,4)<ptrhalt(end,1,4)
    ptrhalt=flip(ptrhalt,1);
end

%calculate mean of the variables for each layer between the sounding levels
layer_p=(ptrhalt(1:size(ptrhalt,1)-1,:,1)+ptrhalt(2:size(ptrhalt,1),:,1))./2;
layer_t=(ptrhalt(1:size(ptrhalt,1)-1,:,2)+ptrhalt(2:size(ptrhalt,1),:,2))./2;
layer_rh=(ptrhalt(1:size(ptrhalt,1)-1,:,3)+ptrhalt(2:size(ptrhalt,1),:,3))./2;

%Find right layers
layer_depth=nan(size(layer_p));
for ii=1:size(ptrhalt,1)-1
    layerUp=ptrhalt(ii,:,4);
    layerDown=ptrhalt(ii+1,:,4);
    
    borderUp=nan(size(layerUp));
    borderDown=nan(size(layerUp));
    
    %in between
    indBetween=find(layerUp>inalt' & layerDown<inalt');
    borderUp(indBetween)=inalt(indBetween);
    borderDown(indBetween)=layerDown(indBetween);
    
    %below flight
    indBelow=find(layerUp<inalt');
    borderUp(indBelow)=layerUp(indBelow);
    borderDown(indBelow)=layerDown(indBelow);
    
    layer_depth(ii,:)=borderUp-borderDown;
end

%saturation pressure after WMO
%layer_es=6.112.*exp((17.62.*layer_t)./(layer_t+243.12));
%saturation pressure after Liebe 1981
layer_es=2.409.*(300./(layer_t+273.15)).^5.*10.^(10-9.834.*(300./(layer_t+273.15))).*10;

%run Liebe code which gives attenuation per kilometer
[spAlpha]=f_atten_Liebe(f,layer_p,layer_t,layer_es,layer_rh);
%multiply by layer depth
layer_alpha=spAlpha.*layer_depth./1000;
%sum over all layers
alphaTot=nansum(layer_alpha,1)';

%run ITU-R code which gives total, dry, and wet, attenuation per kilometer
[spGamma, spGamma0, spGammaW]=f_atten_ITUR(f,layer_p,layer_t,layer_rh,layer_es);
%multiply by layer depth
layer_gamma=spGamma.*layer_depth./1000;
%sum over all layers
gammaTot=nansum(layer_gamma,1)';

end



