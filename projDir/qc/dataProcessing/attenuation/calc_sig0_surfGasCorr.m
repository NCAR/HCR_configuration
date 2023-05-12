function sig0measured=calc_sig0_surfGasCorr(data)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
c=3.0e8;
lambda=c/data.frq;
w_dielec_sq=0.69;

temp=c*pi^5*data.pulse_width*w_dielec_sq/(2*lambda^4*1e18);
rc_terms=10.*log10(temp);

%calculate sig0 measured without bias for all data, not just around
%10 deg
elevTemp=abs(data.elevation+90);
elevTemp(data.elevation>90)=nan;

sig0measured=data.surfRefl-10*log10(cosd(elevTemp))+rc_terms;
end

