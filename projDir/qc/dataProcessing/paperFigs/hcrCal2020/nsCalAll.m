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

wi=10;
hi=4;

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
s1=subplot(1,3,1);
hold on;
in1 = get(ax1,'children');
copyobj(in1,s1);
ylim([-4 0]);
xlim([0,30]);
s1.XTickLabels={'0','10','20',''};
xlabel(['Pod temperature (',char(176),'C)']);
ylabel('DBMVC - noise source (dB)')
text(1,0.2,'(a) No correction','fontsize',12,'fontweight','bold')
grid on

s2=subplot(1,3,2);
hold on;
in2 = get(ax2,'children');
copyobj(in2,s2);
ylim([-4 0]);
xlim([0,30]);
s2.XTickLabels={'0','10','20',''};
s2.YTickLabels='';
xlabel(['Pod temperature (',char(176),'C)']);
text(1,0.2,'(b) LNA temperature corrected','fontsize',12,'fontweight','bold')
grid on

s3=subplot(1,3,3);
hold on;
in3 = get(ax3,'children');
copyobj(in3,s3);
ylim([-4 0]);
xlim([0,30]);
s3.YTickLabels='';
xlabel(['Pod temperature (',char(176),'C)']);
text(1,0.2,'(c) LNA & pod temp. corrected','fontsize',12,'fontweight','bold')
grid on

s1.Position=[0.06 0.12 0.3 0.81];
s2.Position=[0.37 0.12 0.3 0.81];
s3.Position=[0.68 0.12 0.3 0.81];

print(fig1, [figdir,'nsCalAll.png'],'-dpng','-r0');