function [PLT] = f_determine_bias_nadir(PLT,attLiebe,attITU,frq)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
c=299792458;
lambda = c/frq;
w_dielec_sq = PLT.kSquaredWater;

temp = c * pi^5 * PLT.pulseWidth * w_dielec_sq /( 2*lambda^4 * 1e18);
rc_terms = 10 .* log10(temp);

theta = 10.0;  % optimum angle for sigma0 determination

if min(isnan(attLiebe))==0
       
    %calculate sig0 measured without bias for all data, not just around
    %10 deg
    elevTemp=PLT.elev;
    elevTemp(PLT.elev>90)=nan;
    
    PLT.sig0measured=PLT.refl + 2.*attLiebe - 10*log10(cosd(elevTemp)) + rc_terms;
    sig0measuredITU=PLT.refl + 2.*attITU - 10*log10(cosd(elevTemp)) + rc_terms;
    
    % Calculate mean and std of sig0 measured around 10 deg
    sig0measuredMeanL=mean(PLT.sig0measured,'omitnan');
    sig0measuredMeanITU=mean(sig0measuredITU,'omitnan');
       
else
    sig0measuredStd=nan;
end

end

