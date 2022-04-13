% Plot hit miss table

clear all
close all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/paperFigs/'];
indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/comparePID_UW_largest_all/';

sizes=load([indir,'sizes.mat']);
pidSize=sizes.pidSize;
pidStd=sizes.pidStd;

%cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
cscale_hcr=[255,0,0; 255,204,204; 249,163,25; 255,240,60; 136,34,185; 255,0,255; 17,170,51; 0,0,255; 0,255,255; 0,0,0; 150,150,150];

cscale_hcr=cscale_hcr./255;

units_str_hcr={'Rain','SC Rain','Drizzle','SC Drizzle','Cloud Liquid','SC Cloud Liquid',...
    'Melting','Large Frozen','Small Frozen','Precip','Cloud'};

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

errorbar(1:11,pidSize,pidStd,'sk','MarkerSize',0.5,'linewidth',1.5,'capsize',15);
scatter(1:11,pidSize,150,cscale_hcr,'filled','MarkerEdgeColor','k','linewidth',2);
xlim([0 12]);
xticks(1:11);
xticklabels(units_str_hcr);
yticks(0:0.1:2);
ylim([0 1.2]);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel('Mean particle diameter (mm)');

s1.Position=[0.11 0.16 0.86 0.81];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'partSizes.png'],'-dpng','-r0')