%

peak_power_dBm =59.99;
%pulse_width=0.256e-6;
peak_power=10^((peak_power_dBm)*0.1)

wave_len=0.003178e-3; %transmit wavelength in km;

%mod_k=0.711;

tx_ant_gain=10^4.67;
rx_ant_gain=10.^4.67;

%c=300000; % km/sec
az_beam=0.75;
el_beam=0.72; % in degrees

rc1=peak_power*(tx_ant_gain*rx_ant_gain)*wave_len*wave_len;
rc3=az_beam*el_beam*pi*pi/(180*180);
rc4= 512*log(2)*pi*pi;
Wband_RC=10*log10(rc1*rc3/rc4)

Pr=-50;%dBm

Range=0.5:0.5:15; % in km
loss=1.7; % Rx mismatch loss=1.1 dB and Two-way radome loss=0.6 dB
avg_WV_loss= 3; % Two-way WV attenuation in  dB
%sea_surface_bias=4.6+2.1; %CSET=4.6 dB; SOCRATES=2.1 dB; Measured dBZ is lower by 4.6+2.1 dB;
Wband_RC_revised=Wband_RC-loss-avg_WV_loss;
dBZW=-Wband_RC_revised+Pr+20*log10(Range)-10*log10(cosd(10));
plot(Range,dBZW,'b');
%title('Sensitivity of Radar','fontsize',14,'fontweight','b');
    xlabel('Range, km','fontsize',16,'fontweight','b');
    ylabel('Sigma0, dB','fontsize',16,'fontweight','b');

function g_abs=absorb(ele_deg,r_km)
fac1=0.4+3.45*exp(-ele_deg/1.8);
fac2=1-exp(-r_km/(27.8+154*exp(-ele_deg/2.2)));
g_abs=fac1*fac2
end
