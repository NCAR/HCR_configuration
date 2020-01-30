% Calibration paper figure of theoretical sig0 models

clear all;
close all;

savefig=1;

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');
addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

figdir='/h/eol/romatsch/hcrCalib/oceanScans/figsCalibPaper/';

PLT.sig0measured=0;
PLT.elev=0:0.1:20;
PLT.elev=PLT.elev';

%% Vary wind
oceanTemp=20;
salinity=35;
frq=9.4406e+10;

surface_wind1={2};
[PLT1]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

surface_wind1={5};
[PLT2]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

surface_wind1={10};
[PLT3]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

surface_wind1={15};
[PLT4]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

surface_wind1={20};
[PLT5]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

%% Vary SST
surface_wind1={5};
salinity=35;
frq=9.4406e+10;

oceanTemp=0;
[PLT6]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

oceanTemp=10;
[PLT7]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

oceanTemp=20;
[PLT8]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);

oceanTemp=30;
[PLT9]= f_sigma0_model(PLT,surface_wind1,frq,oceanTemp,salinity);


%% Figure
close all

wi=5;
hi=8;

fig=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[1,100,wi,hi]);
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPosition = [0, 0, wi, hi];
fig.PaperSize = [wi, hi];
fig.Resize = 'off';
fig.InvertHardcopy = 'off';

set(fig,'color','w');

ax1=subplot(2,1,1,'units','inch');
hold on;
ax1.Position = [0.6 4.55 4.0 3.1];

cols=jet(5);

plot(PLT.elev,PLT1.sig0model(:,8),'color',cols(1,:),'linewidth',2);
plot(PLT.elev,PLT2.sig0model(:,8),'color',cols(2,:),'linewidth',2);
plot(PLT.elev,PLT3.sig0model(:,8),'color',cols(3,:),'linewidth',2);
plot(PLT.elev,PLT4.sig0model(:,8),'color',cols(4,:),'linewidth',2);
plot(PLT.elev,PLT5.sig0model(:,8),'color',cols(5,:),'linewidth',2);

xlim([0,20]);
ylim([-5 15]);
grid on
xlabel('Pointing angle off nadir [deg]');
ylabel('sig0 [dB]');

legend('2 m s^{-1}','5 m s^{-1}','10 m s^{-1}','15 m s^{-1}','20 m s^{-1}');
title('a           Variation with surface wind speed')

ax2=subplot(2,1,2,'units','inch');
hold on;
ax2.Position = [0.6 0.55 4.0 3.1];

cols=jet(5);

plot(PLT.elev,PLT6.sig0model(:,8),'color',cols(1,:),'linewidth',2);
plot(PLT.elev,PLT7.sig0model(:,8),'color',cols(2,:),'linewidth',2);
plot(PLT.elev,PLT8.sig0model(:,8),'color',cols(3,:),'linewidth',2);
plot(PLT.elev,PLT9.sig0model(:,8),'color',cols(5,:),'linewidth',2);

xlim([0,20]);
ylim([-5 15]);
grid on
xlabel('Pointing angle off nadir [deg]');
ylabel('sig0 [dB]');

legend('0 C','10 C','20 C','30 C');
title('b  Variation with sea surface temperature')

if savefig
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'CoxMunk_windSST'],'-dpng','-r0');
end
