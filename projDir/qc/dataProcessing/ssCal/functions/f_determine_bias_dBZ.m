function [PLT dB_bias_lb dB_bias_ITU sig0measuredStd N] = f_determine_bias_dBZ(PLT,attLiebe,attITU,frq)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
c = 3.0e8;
tau = 2.56e-7;
lambda = c/frq;
w_dielec_sq = .69;

temp = c * pi^5 * tau * w_dielec_sq /( 2*lambda^4 * 1e18);
rc_terms = 10 * log10(temp);

theta = 10.0;  % optimum angle for sigma0 determination
dB_sigma10 = 6.0;   % expected sigma0 at 10 degrees

% find beams near 10-deg elevation angle:
A = find( abs(theta - PLT.elev) < 0.3);

if ~isnan(attLiebe)
       
    %calculate sig0 measured without bias for all data, not just around
    %10 deg
    PLT.sig0measured=PLT.refl + 2*attLiebe - 10*log10(cosd(PLT.elev)) + rc_terms;
    sig0measuredITU=PLT.refl + 2*attITU - 10*log10(cosd(PLT.elev)) + rc_terms;
    
    % Calculate mean and std of sig0 measured around 10 deg
    sig0measuredMeanL=nanmean(PLT.sig0measured(A));
    sig0measuredMeanITU=nanmean(sig0measuredITU(A));
    sig0measuredStd=nanstd(PLT.sig0measured(A));
    N=length(A);
    
    % Calculate mean bias
    dB_bias_lb=dB_sigma10-sig0measuredMeanL;   
    dB_bias_ITU=dB_sigma10-sig0measuredMeanITU;
else
    dB_bias_lb=nan;
    dB_bias_ITU=nan;
    PLT.sig0measured=nan;
    sig0measuredStd=nan;
    N=nan;
end

end

