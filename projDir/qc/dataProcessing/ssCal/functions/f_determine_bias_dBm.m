function [PLT dB_bias_lb dB_bias_ITU sig0measuredStd N] = f_determine_bias_dBm(PLT,attLiebe,attITU,frq)

peak_power_dBm =PLT.xmitPowV;
peak_power=10^((peak_power_dBm)*0.1);

c = 3.0e8;
wave_len = c/frq/1000;%transmit wavelength in km;

tx_ant_gain=10^(PLT.antGainV/10);

az_beam=PLT.beamWidthH;
el_beam=PLT.beamWidthV; % in degrees

rc1=peak_power*(tx_ant_gain*tx_ant_gain)*wave_len*wave_len;
rc3=az_beam*el_beam*pi*pi/(180*180);
rc4= 512*log(2)*pi*pi;
Wband_RC=10*log10(rc1*rc3/rc4);

Pr=PLT.vpwr;%dBm
Range=PLT.range./1000; % in km

loss=PLT.waveGuideLoss+PLT.radomeLoss+PLT.recMismatchLoss;

theta = 10.0;  % optimum angle for sigma0 determination
dB_sigma10 = 6.0;   % expected sigma0 at 10 degrees

% find beams near 10-deg elevation angle:
A = find( abs(theta - PLT.elev) < 0.3);

if ~isnan(attLiebe)
    
    %calculate sig0 measured without bias for all data, not just around
    %10 deg
    avg_WV_lossL= 2*attLiebe; % Two-way WV attenuation in  dB
    Wband_RC_revisedL=Wband_RC-loss-avg_WV_lossL;
    PLT.sig0measured=-Wband_RC_revisedL+Pr+20*log10(Range)-10*log10(cosd(PLT.elev));
    
    avg_WV_lossITU= 2*attITU; % Two-way WV attenuation in  dB
    Wband_RC_revisedITU=Wband_RC-loss-avg_WV_lossITU;
    sig0measuredITU=-Wband_RC_revisedITU+Pr+20*log10(Range)-10*log10(cosd(PLT.elev));
    
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

% plot(PLT.elev,PLT.sig0measured,'b');
% title(datestr(PLT.time(1)),'fontsize',14,'fontweight','b');
% xlabel('Elev, deg','fontsize',16,'fontweight','b');
% ylabel('Sigma0, dB','fontsize',16,'fontweight','b');
end
