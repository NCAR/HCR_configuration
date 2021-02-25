function sig0model= calc_sig0_model(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lihua Li et al 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
salinity=35;
surfWind=sqrt(data.U_SURF.^2+data.V_SURF.^2);

if max(~isnan(surfWind))~=0
        
    sig0model=nan(9,size(data.elevation,2));
    
    c = 3.0e8;
    lambda = c/data.frq;
    
    n=f_refractiveIndex_MeissnerWentz(data.SST,salinity,data.frq);
    
    Ce=0.88; %Li et al. 2005
    CeLow=0.85;
    CeHigh=0.95;
    
    GammaE=Ce.*(n-1)./(n+1);
    GammaELow=CeLow.*(n-1)./(n+1);
    GammaEHigh=CeHigh.*(n-1)./(n+1);
    
    % Freilich and Vanhoff 2003
    w0FV=nan(size(data.elevation));
    w1FV=nan(size(data.elevation));
    
    w0FV(surfWind>1 & surfWind<=10)=0.0036;
    w1FV(surfWind>1 & surfWind<=10)=0.028;
    w0FV(surfWind>10 & surfWind<20)=-0.0184;
    w1FV(surfWind>10 & surfWind<20)=0.05;
    
    s2FV=w0FV+w1FV.*log10(surfWind);
    
    sig0model(1,:)=10.*log10(abs(GammaELow).^2./(s2FV.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2FV)));
    sig0model(2,:)=10.*log10(abs(GammaE).^2./(s2FV.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2FV)));
    sig0model(3,:)=10.*log10(abs(GammaEHigh).^2./(s2FV.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2FV)));
    
    % Wu 1972
    w0Wu=nan(size(data.elevation));
    w1Wu=nan(size(data.elevation));
    
    w0Wu(surfWind>1 & surfWind<=7)=0.009;
    w1Wu(surfWind>1 & surfWind<=7)=0.0276;
    w0Wu(surfWind>7 & surfWind<20)=-0.084;
    w1Wu(surfWind>7 & surfWind<20)=0.138;
    
    s2Wu=w0Wu+w1Wu.*log10(surfWind);
    
    sig0model(4,:)=10.*log10(abs(GammaELow).^2./(s2Wu.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2Wu)));
    sig0model(5,:)=10.*log10(abs(GammaE).^2./(s2Wu.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2Wu)));
    sig0model(6,:)=10.*log10(abs(GammaEHigh).^2./(s2Wu.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2Wu)));
    
    % Cox and Munk 1954
    s2CM=0.003+5.08e-3.*surfWind;
    
    sig0model(7,:)=10.*log10(abs(GammaELow).^2./(s2CM.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2CM)));
    sig0model(8,:)=10.*log10(abs(GammaE).^2./(s2CM.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2CM)));
    sig0model(9,:)=10.*log10(abs(GammaEHigh).^2./(s2CM.*(cosd(data.elevation)).^4).*exp(-(tand(data.elevation)).^2./(s2CM)));
    
else
    sig0model=nan;
end

end



