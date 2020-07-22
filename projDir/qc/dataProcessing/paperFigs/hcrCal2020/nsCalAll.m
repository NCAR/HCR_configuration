% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

%% Load figures
h1 = openfig([figdir,'nscal/lnaTemps.fig'],'reuse');
ax1 = gca;
h2 = openfig([figdir,'nscal/podTemps_line.fig'],'reuse');
ax2 = gca;
h3 = openfig([figdir,'nscal/checkTemp.fig'],'reuse');
pause(1);
ax3 = gca;

%% Plot

wi=5;
hi=12;

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
s1=subplot(3,1,1);
hold on;
outerpos1 = s1.Position;
s1.Position = [outerpos1(1)-0.0 outerpos1(2)+0.0 outerpos1(3)+0.03 outerpos1(4)+0.05];
in1 = get(ax1,'children');
copyobj(in1,s1);
ylim([-4 0]);
xlim([0,30]);
xlabel('Pod temperature (C)');
ylabel('DBMVC - noise source (dB)')
text(1,0.2,'(a) No correction','fontsize',12,'fontweight','bold')
grid on

s2=subplot(3,1,2);
hold on;
outerpos1 = s2.Position;
s2.Position = [outerpos1(1)-0.0 outerpos1(2)-0.035 outerpos1(3)+0.03 outerpos1(4)+0.05];
in2 = get(ax2,'children');
copyobj(in2,s2);
ylim([-4 0]);
xlim([0,30]);
xlabel('Pod temperature (C)');
ylabel('DBMVC - noise source (dB)')
text(1,0.2,'(b) LNA temperature corrected','fontsize',12,'fontweight','bold')
grid on

s3=subplot(3,1,3);
hold on;
outerpos1 = s3.Position;
s3.Position = [outerpos1(1)-0.0 outerpos1(2)-0.07 outerpos1(3)+0.03 outerpos1(4)+0.05];
in3 = get(ax3,'children');
copyobj(in3,s3);
ylim([-4 0]);
xlim([0,30]);
xlabel('Pod temperature (C)');
ylabel('DBMVC - noise source (dB)')
text(1,0.2,'(c) LNA and pod temperature corrected','fontsize',12,'fontweight','bold')
grid on

print(fig1, [figdir,'nsCalAll.png'],'-dpng','-r0');