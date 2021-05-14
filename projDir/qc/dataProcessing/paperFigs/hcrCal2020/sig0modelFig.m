% Calibration paper figure of theoretical sig0 models

clear all;
close all;

savefig=1;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

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

wi=10;
hi=5;

fig=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[1,100,wi,hi],'renderer','painters');
fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'inches';
fig.Units = 'inches';
fig.PaperPosition = [0, 0, wi, hi];
fig.PaperSize = [wi, hi];
fig.Resize = 'off';
fig.InvertHardcopy = 'off';

set(fig,'color','w');

ax1=subplot(1,2,1,'units','inch');
hold on;

cols=jet(5);

plot(PLT.elev,PLT1.sig0model(:,8),'color',cols(1,:),'linewidth',2);
plot(PLT.elev,PLT2.sig0model(:,8),'color',cols(2,:),'linewidth',2);
plot(PLT.elev,PLT3.sig0model(:,8),'color',cols(3,:),'linewidth',2);
plot(PLT.elev,PLT4.sig0model(:,8),'color',cols(4,:),'linewidth',2);
plot(PLT.elev,PLT5.sig0model(:,8),'color',cols(5,:),'linewidth',2);

xlim([0,20]);
ylim([-5 15]);
grid on
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)');

legend('2 m s^{-1}','5 m s^{-1}','10 m s^{-1}','15 m s^{-1}','20 m s^{-1}');
text(0,16,'(a) \sigma_0 variation with surface wind speed','fontsize',12,'fontweight','bold')

ax2=subplot(1,2,2,'units','inch');
hold on;

cols=jet(5);

plot(PLT.elev,PLT6.sig0model(:,8),'color',cols(1,:),'linewidth',2);
plot(PLT.elev,PLT7.sig0model(:,8),'color',cols(2,:),'linewidth',2);
plot(PLT.elev,PLT8.sig0model(:,8),'color',cols(3,:),'linewidth',2);
plot(PLT.elev,PLT9.sig0model(:,8),'color',cols(5,:),'linewidth',2);

xlim([0,20]);
ylim([-5 15]);
grid on
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)');

legend(['0 ',char(176),'C'],['10 ',char(176),'C'],['20 ',char(176),'C'],['30 ',char(176),'C']);
text(0,16,'(b) \sigma_0 variation with SST','fontsize',12,'fontweight','bold')

ax1.Position = [0.65 0.55 4.2 4];
ax2.Position = [5.6 0.55 4.2 4];

if savefig
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'CoxMunk_windSST'],'-dpng','-r0');
end
