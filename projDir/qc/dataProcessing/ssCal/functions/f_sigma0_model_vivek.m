function [PLT]= f_processOceanScans_dBm(PLT,surface_wind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lihua Li Seattle AMS Radar conference  pp 204 -207, 6-12, August 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if min(isnan(PLT.sig0measured))==0 & ~isnan(surface_wind)
    
    ref_index= complex(7.56, -13.6);  % Refractive index for 20 C sea surafce temperature
    
    eff_ref_coeff= 0.89*(ref_index-1)/(ref_index+1);  % 0.89 Ocean surface roughness correction factor Freilich, 2003
    
    fac1=abs(eff_ref_coeff)^2;
        
    elangle=PLT.elev*pi/180;
    sur_slope= (0.003+1.92e-3*surface_wind)+ (3.16e-3*surface_wind); % Cox and Munk 1954
    
    fac2= sur_slope.*(cos(elangle)).^4;
    fac3=exp(-(tan(elangle)).^2./sur_slope);
        
    PLT.sig0model=10.*log10(((fac1*fac3./fac2)));
else
    PLT.sig0model=nan;
end

end



