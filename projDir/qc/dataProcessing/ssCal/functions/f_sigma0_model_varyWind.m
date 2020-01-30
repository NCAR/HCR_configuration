function [PLT]= f_sigma0_model_varyWind(PLT,surface_wind,frq,oceanTemp,salinity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lihua Li et al 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varWind=5;
surfWind=surface_wind{:};
surfWindLow=surfWind-varWind;
surfWindHigh=surfWind+varWind;

surfWindAll=cat(2,surfWindLow,surfWind,surfWindHigh);

PLT.sig0model=nan(size(PLT.sig0measured,1),3);

for ii=1:size(surfWindAll,2)
    surfWind=surfWindAll(:,ii)
    
    if min(isnan(PLT.sig0measured))==0 & max(~isnan(surfWind))~=0
        
        if size(surfWind,1)~=size(PLT.time,1)
            surfWindTemp=nan(size(PLT.time));
            surfWindTemp(:)=surfWind;
            surfWind=surfWindTemp;
        end
        
        c = 3.0e8;
        lambda = c/frq;
        
        n=f_refractiveIndex_MeissnerWentz(oceanTemp,salinity,frq);
        
        Ce=0.88; %Li et al. 2005
        CeLow=0.85;
        CeHigh=0.95;
        
        GammaE=Ce.*(n-1)./(n+1);
                
        % Cox and Munk 1954
        s2CM=0.003+5.08e-3.*surfWind;
        
        PLT.sig0model(:,ii)=10.*log10(abs(GammaE).^2./(s2CM.*(cosd(PLT.elev)).^4).*exp(-(tand(PLT.elev)).^2./(s2CM)));
               
    else
        PLT.sig0model=nan;
    end
end
end



