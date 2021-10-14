% Plot hit miss table

clear all
close all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/paperFigs/'];
indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/comparePID_UW_wholeFlights/';

sizes=load([indir,'sizes.mat']);
pidSize=sizes.pidSize;
pidStd=sizes.pidStd;

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
units_str_hcr={'Rain','SC rain','Drizzle','SC drizzle','Cloud liquid','SC cloud liquid','Mixed phase','Large frozen','Small frozen'};

wi=5.2;
hi=5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi],'renderer','painters');
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(1,1,1);
hold on

errorbar(1:9,pidSize,pidStd,'sk','MarkerSize',0.5,'linewidth',1.5,'capsize',15);
scatter(1:9,pidSize,150,cscale_hcr,'filled','MarkerEdgeColor','k','linewidth',2);
xlim([0 10]);
xticks(1:9);
xticklabels(units_str_hcr);
yticks(0:0.1:2);
ylim([0 1.6]);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel('Mean particle size (mm)');

s1.Position=[0.11 0.16 0.86 0.81];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'partSizes.png'],'-dpng','-r0')