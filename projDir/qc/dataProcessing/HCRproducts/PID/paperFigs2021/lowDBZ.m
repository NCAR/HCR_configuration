% Plot hit miss table

clear all
close all

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/paperFigs/'];

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/lowRefl/';

outTableIn=load([indir,'pid_dbzLow.mat']);
dbzAll=outTableIn.dbzAll;
pidAll=outTableIn.pidAll;

%cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
cscale_hcr=[255,0,0; 255,204,204; 249,163,25; 255,240,60; 136,34,85; 255,0,255; 17,170,51; 0,0,255; 0,255,255];

cscale_hcr=cscale_hcr./255;

dbzPID=cat(2,dbzAll,pidAll);
dbzPID(find(dbzAll>-30),:)=[];

edges=0.5:9.5;
countsCat=histcounts(dbzPID,edges);

units_str_hcr={'Rain','SC Rain','Drizzle','SC Drizzle','Cloud Liquid','SC Cloud Liquid',...
    'Melting','Large Frozen','Small Frozen','Precip','Cloud'};

wi=6;
hi=4.6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(1,1,1);

b=bar(1:9,countsCat./sum(countsCat).*100,1);

b.FaceColor='flat';
for jj=1:9
    b.CData(jj,:)=cscale_hcr(jj,:);
end

xlim([0.5,9.5]);
ylim([0 70]);

xticklabels(units_str_hcr);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
ylabel('Particles (%)')

sumLiq=sum(countsCat(1:6));
sumFrozen=sum(countsCat(8:9));

cloudPerc=sumLiq./sum(countsCat).*100;
sFr=sumFrozen./sum(countsCat).*100;

text(1,66,['Liquid: ',num2str(cloudPerc,2),' %. Frozen: ',num2str(sFr,2),' %.'],'fontsize',12,'FontWeight','bold');

s1.Position=[0.09 0.18 0.9 0.8];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'lowRefl.png'],'-dpng','-r0')
