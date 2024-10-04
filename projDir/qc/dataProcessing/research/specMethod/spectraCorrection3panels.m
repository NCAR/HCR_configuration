clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'spectraCorr2.mat']);

%% Figure

close all
xlim1=[1,988];
xlim2=[460,562];

f1 = figure('Position',[200 500 600 1000],'DefaultAxesFontSize',12,'renderer','painters');

t = tiledlayout(3,1,'TileSpacing','tight','Padding','compact','TileIndexing', 'columnmajor');

s1=nexttile(1);

hold on

l1=plot(spectra.rawDB,'-b','LineWidth',1);
l2=plot(xlim1,[spectra.noiseFloor,spectra.noiseFloor],'-c','LineWidth',1.5);
ylabel('Power (dB)');
ylim([-80,-20]);
xlim(xlim1);

yyaxis right
l3=plot(spectra.BBS,'-r','LineWidth',2);
xlabel('Spectral bin number');
set(gca,'YColor','r');
ylim([0,1.2]);

grid on
box on
%title('(a) Example spectrum')
legend([l1,l3,l2],{'PS_{raw}';'GAUSS_{bbs}';'Spec. noise floor'},'Location','northwest');
text(920,1.1,'(a)','FontSize',12,'FontWeight','bold')

s2=nexttile(2);
dftps=cat(2,spectra.fftRaw(495:end),spectra.fftRaw(1:494));
dftbbs=cat(2,spectra.fftBBS(495:end),spectra.fftBBS(1:494));

xdata=19:1006;

hold on
l4=plot(xdata,dftbbs,'-r','LineWidth',2);
s2.YAxisLocation = 'right';
s2.XTick = [];
s2.YAxis.Color = 'r';
xlim(xlim2);
yticks(0:0.025:0.1);

s22=axes(t,'Color','none');
s22.Layout.Tile=2;
hold on
l5=plot(xdata,dftps,'-b','LineWidth',2);
xlim(xlim2);
xlabel('Spectral bin number');
ylim([0,60]);
yticks(0:15:75);
annotation('textarrow',[0.65,0.56],[0.45,0.38],'String','Truncation point','FontSize',12)

grid on
box on
%title('(a) Example spectrum')
legend([l5,l4],{'DFT_{ps}';'DFT_{bbs}'},'Location','northwest');
text(555,55,'(b)','FontSize',12,'FontWeight','bold')

s3=nexttile(3);

hold on

l1=plot(spectra.rawDB,'-b','LineWidth',1);
l2=plot(xlim1,[spectra.noiseFloor,spectra.noiseFloor],'-c','LineWidth',2);
l3=plot(spectra.uncorrFilteredDB,'-','LineWidth',2,'Color',[0.8,0,0.8]);
l4=plot(spectra.corrDB,'-g','LineWidth',2);
ylabel('Power (dB)');
xlabel('Spectral bin number');
ylim([-80,-20]);
xlim(xlim1);

grid on
box on
%title('(a) Example spectrum')
legend([l1,l4,l3,l2],{'PS_{raw}';'PS_{cf}';'PS_{f}';'Sp. noise floor'},'Location','northwest');
text(920,-25,'(c)','FontSize',12,'FontWeight','bold')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'specCorr.png'],'-dpng','-r0');
