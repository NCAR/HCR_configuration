% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];
figdirIn=['/h/eol/romatsch/hcrCalib/sensitivity/'];

%% Load figures
h1 = openfig([figdirIn,'otrec_minDBZ_1km.fig'],'reuse');
in1=subplot(2,2,2);
ax1 = gca;

in2=subplot(2,2,4);
ax2 = gca;

%% Plot

wi=5;
hi=9;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1=subplot(2,1,1);
hold on;
outerpos1 = s1.Position;
s1.Position = [outerpos1(1)-0.0 outerpos1(2)+0.0 outerpos1(3)+0.03 outerpos1(4)+0.05];
inget1 = get(ax1,'children');
copyobj(inget1,s1);
%ylim([-4 0]);
xlim([-50,-35]);
ylabel('Data count');
xlabel('Reflectivity (dBZ)')
text(-49,3.8e+4,'(a) DBZ','fontsize',12,'fontweight','bold')
grid on

s2=subplot(2,1,2);
hold on;
outerpos1 = s2.Position;
s2.Position = [outerpos1(1)-0.0 outerpos1(2)-0.035 outerpos1(3)+0.03 outerpos1(4)+0.05];
inget2 = get(ax2,'children');
copyobj(inget2,s2);
%ylim([-4 0]);
xlim([-20,-5]);
ylabel('Data count');
xlabel('Signal to Noise Ratio (dB)')
text(-19,5.7e+4,'(b) SNR','fontsize',12,'fontweight','bold')
grid on

print(fig1, [figdir,'sensitivity.png'],'-dpng','-r0');