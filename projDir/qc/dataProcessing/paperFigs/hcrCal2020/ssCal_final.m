% Correct HCR altitude
clear all;
close all;

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

%% Load figures
h1 = openfig([figdir,'cset_qc2_file_sig0model_paper/cset_qc2_file_sig0model_paper_bias1.fig'],'reuse');
ax1 = gca;
h2 = openfig([figdir,'socrates_qc2_era5file_sig0model_paper/socrates_qc2_era5file_sig0model_paper_bias1.fig'],'reuse');
ax2 = gca;
h3 = openfig([figdir,'otrec_qc2_era5file_sig0model_paper/otrec_qc2_era5file_sig0model_paper_bias1.fig'],'reuse');
ax3 = gca;
h4 = openfig([figdir,'cset_qc2_file_sig0model_paper/cset_qc2_file_sig0model_paper_bias2.fig'],'reuse');
ax4 = gca;
h5 = openfig([figdir,'socrates_qc2_era5file_sig0model_paper/socrates_qc2_era5file_sig0model_paper_bias2.fig'],'reuse');
ax5 = gca;
h6 = openfig([figdir,'otrec_qc2_era5file_sig0model_paper/otrec_qc2_era5file_sig0model_paper_bias2.fig'],'reuse');
pause(1)
ax6 = gca;

%% Plot
%close
wi=10;
hi=6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[690,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1=subplot(2,3,1);
hold on;
outerpos1 = s1.Position;
s1.Position = [outerpos1(1)-0.075 outerpos1(2)-0.0 0.3 0.35];
in1 = get(ax1,'children');
copyobj(in1,s1);
ylim([-4 14]);
xlim([5,15]);
xticks(6:15)
yticks(-4:2:14);
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)')
text(5,15.5,'(a) CSET','fontsize',12,'fontweight','bold')
grid on

s2=subplot(2,3,2);
hold on;
outerpos1 = s2.Position;
s2.Position = [outerpos1(1)-0.04 outerpos1(2)-0.0 0.3 0.35];
in2 = get(ax2,'children');
copyobj(in2,s2);
ylim([-4 14]);
xlim([5,15]);
xticks(6:15);
yticks(-4:2:14);
yticklabels('');
xlabel('Incidence angle (deg)');
text(5,15.5,'(b) SOCRATES','fontsize',12,'fontweight','bold')
grid on

s3=subplot(2,3,3);
hold on;
outerpos1 = s3.Position;
s3.Position = [outerpos1(1)-0.005 outerpos1(2)-0.0 0.3 0.35];
in3 = get(ax3,'children');
copyobj(in3,s3);
ylim([-4 14]);
xlim([5,15]);
xticks(6:15)
yticks(-4:2:14);
yticklabels('');
xlabel('Incidence angle (deg)');
text(5,15.5,'(c) OTREC','fontsize',12,'fontweight','bold')
grid on

s4=subplot(2,3,4);
hold on;
outerpos1 = s4.Position;
s4.Position = [outerpos1(1)-0.075 outerpos1(2)-0.02 0.3 0.38];
in4 = get(ax4,'children');
plot([5,15],[0,0],'-k')
copyobj(in4,s4);
ylim([-4 7]);
xlim([5,15]);
xticks(6:15)
yticks(-4:1:14);
xlabel('Incidence angle (deg)');
ylabel('\sigma_0 (dB)')
text1= findobj(gcf, 'Type', 'Text');
line1=text1(1);
line1.FontSize=11;
line1.String = strrep(line1.String,' model','');
line1.String = strrep(line1.String,' bias','');
line2=text1(2);
line2.FontSize=11;
line2.String = strrep(line2.String,' model','');
line2.String = strrep(line2.String,' bias','');
line3=text1(3);
line3.FontSize=11;
line3.String = strrep(line3.String,' model','');
line3.String = strrep(line3.String,' bias','');
text(5,8,'(d)','fontsize',12,'fontweight','bold')
grid on

s5=subplot(2,3,5);
hold on;
outerpos1 = s5.Position;
s5.Position = [outerpos1(1)-0.04 outerpos1(2)-0.02 0.3 0.38];
in5 = get(ax5,'children');
plot([5,15],[0,0],'-k')
copyobj(in5,s5);
ylim([-4 7]);
xlim([5,15]);
xticks(6:15);
yticks(-4:1:14);
yticklabels('');
xlabel('Incidence angle (deg)');
text1= findobj(gcf, 'Type', 'Text');
line1=text1(1);
line1.FontSize=11;
line1.String = strrep(line1.String,' model','');
line1.String = strrep(line1.String,' bias','');
line2=text1(2);
line2.FontSize=11;
line2.String = strrep(line2.String,' model','');
line2.String = strrep(line2.String,' bias','');
line3=text1(3);
line3.FontSize=11;
line3.String = strrep(line3.String,' model','');
line3.String = strrep(line3.String,' bias','');
text(5,8,'(e)','fontsize',12,'fontweight','bold')
grid on

s6=subplot(2,3,6);
hold on;
outerpos1 = s6.Position;
s6.Position = [outerpos1(1)-0.005 outerpos1(2)-0.02 0.3 0.38];
in6 = get(ax6,'children');
plot([5,15],[0,0],'-k')
copyobj(in6,s6);
ylim([-4 7]);
xlim([5,15]);
xticks(6:15)
yticks(-4:1:14);
yticklabels('');
xlabel('Incidence angle (deg)');
text1= findobj(gcf, 'Type', 'Text');
line1=text1(1);
line1.FontSize=11;
line1.String = strrep(line1.String,' model','');
line1.String = strrep(line1.String,' bias','');
line2=text1(2);
line2.FontSize=11;
line2.String = strrep(line2.String,' model','');
line2.String = strrep(line2.String,' bias','');
line3=text1(3);
line3.FontSize=11;
line3.String = strrep(line3.String,' model','');
line3.String = strrep(line3.String,' bias','');
text(5,8,'(f)','fontsize',12,'fontweight','bold')
grid on

allLines=findall(gcf, 'type', 'line');
l1=allLines(3);
l2=allLines(2);
l3=allLines(1);
leg=legend([l1,l2,l3],{'FV model','Wu model','CM model'},'location','southwest','fontsize',11);
leg.Position=[0.69    0.09    0.1    0.1024];
leg.ItemTokenSize = [13,18];
set(gcf,'children',flipud(get(gcf,'children')))

print(fig1, [figdir,'ssCal_final.png'],'-dpng','-r0');