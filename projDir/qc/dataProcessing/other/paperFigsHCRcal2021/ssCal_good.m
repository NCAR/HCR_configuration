% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

%% Load figures
h1 = openfig([figdir,'cset_qc2_file_sig0model_paper/sig0measured/cset_qc2_sig0model_paper_sig0measured_20150709_222330.fig'],'reuse');
ax1 = gca;
h2 = openfig([figdir,'otrec_qc2_era5file_sig0model_paper/sig0measured/otrec_qc2_sig0model_paper_sig0measured_20190818_133600.fig'],'reuse');
pause(1);
ax2 = gca;

%% Plot

wi=10;
hi=5;

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
s1=subplot(1,2,1);
hold on;
in1 = get(ax1,'children');
copyobj(in1,s1);
ylim([-20 15]);
xlim([0,25]);
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)')
text1= findobj(gcf, 'Type', 'Text');
firstLine=text1(1);
firstLine.Position=[10 12 0];
secondLine=text1(2);
secondLine.Position=[10 14 0];
text(0.1,17,'(a) CSET 2015-07-09 22:23:30','fontsize',12,'fontweight','bold')
grid on

s2=subplot(1,2,2);
hold on;
in2 = get(ax2,'children');
copyobj(in2,s2);
ylim([-20 15]);
xlim([0,25]);
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)')
text1= findobj(gcf, 'Type', 'Text');
firstLine=text1(1);
firstLine.Position=[10 12 0];
secondLine=text1(2);
secondLine.Position=[10 14 0];
text(0.1,17,'(b) OTREC 2019-08-18 13:36:00','fontsize',12,'fontweight','bold')
grid on
allLines=findall(gcf, 'type', 'line');
l1=allLines(14);
l2=allLines(13);
l3=allLines(12);
l4=allLines(11);
l5=allLines(10);
l7=allLines(8);
leg=legend([l1,l2,l3,l4,l5,l7],{'Rot angle <= 180 deg','Rot angle > 180 deg','FV model','Wu model','CM model','Data fit'},'location','southwest','fontsize',12);
set(gcf,'children',flipud(get(gcf,'children')))

s1.Position=[0.07 0.11 0.42 0.81];
s2.Position=[0.56 0.11 0.42 0.81];

leg.Position=[0.58 0.2 0.27 0.15];

print(fig1, [figdir,'ssCal_good.png'],'-dpng','-r0');