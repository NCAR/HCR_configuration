% Plot hit miss table

clear all
close all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/paperFigs/'];
indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/comparePID_UW_wholeFlights/';

hmTableIn=load([indir,'hitMissTable.mat']);
hmTableL=hmTableIn.hmTableL;

hmNormL=hmTableL./sum(sum(hmTableL)).*100;

xvalues = {'Frozen','Mixed','Liquid'};

wi=5.2;
hi=5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

h=heatmap(xvalues,xvalues,hmNormL);
ax = gca;
axp = struct(ax);       %you will get a warning
axp.Axes.XAxisLocation = 'top';
h.ColorbarVisible = 'off';

h.XLabel = 'PID';
h.YLabel = 'UWILD';
h.CellLabelFormat = '%.1f';
h.FontSize=11;
h.Position=[0.17 0.06 0.8 0.84];
axes('visible','off')
text(0.1,-0.1,['Correct: ',num2str(hmNormL(1,1)+hmNormL(2,2)+hmNormL(3,3),3),'%. Correlation coefficient: 0.65.'],...
    'fontsize',11);

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'hitMiss.png'],'-dpng','-r0')