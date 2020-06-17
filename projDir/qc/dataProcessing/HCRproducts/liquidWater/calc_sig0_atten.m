function [data] = calc_sig0_atten(data,frq)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
c = 3.0e8;
lambda = c/frq;
w_dielec_sq = .69;

temp = c * pi^5 * data.pulseWidth * w_dielec_sq /( 2*lambda^4 * 1e18);
rc_terms = 10 .* log10(temp);

theta = 10.0;  % optimum angle for sigma0 determination

%calculate sig0 measured without bias for all data, not just around
%10 deg
elevTemp=data.elev;
elevTemp(data.elev>90)=nan;

data.sig0measured=data.refl - 10*log10(cosd(elevTemp)) + rc_terms;
end

