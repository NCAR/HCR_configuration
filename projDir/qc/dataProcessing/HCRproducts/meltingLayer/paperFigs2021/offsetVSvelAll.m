% Analyze HCR clouds

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
indir=['/home/romatsch/plots/HCR/meltingLayer/offsetVel/'];

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

%% Load data
cset=load([indir,'cset_velOffset.mat']);
offsetC=cset.offset;
velAboveC=cset.velAbove;
velBelowC=cset.velBelow;

socrates=load([indir,'socrates_velOffset.mat']);
offsetS=socrates.offset;
velAboveS=socrates.velAbove;
velBelowS=socrates.velBelow;

otrec=load([indir,'otrec_velOffset.mat']);
offsetO=otrec.offset;
velAboveO=otrec.velAbove;
velBelowO=otrec.velBelow;

offset=cat(1,offsetC,offsetS,offsetO);
velAbove=cat(1,velAboveC,velAboveS,velAboveO);
velBelow=cat(1,velBelowC,velBelowS,velBelowO);

%% Plot

close all
edges={-10:0.2:10 -100:20:800};

N=hist3(cat(2,velAbove,offset),'Edges',edges);
% Regression
fitOrthA=gmregress(velAbove,offset,1);
fitAllA=[fitOrthA(2) fitOrthA(1)];
xFitA = -10:0.2:10;
yFitA = polyval(fitAllA, xFitA);

N2=hist3(cat(2,velBelow,offset),'Edges',edges);
% Regression
fitOrthB=gmregress(velBelow,offset,1);
fitAllB=[fitOrthB(2) fitOrthB(1)];
xFitB = -10:0.2:10;
yFitB = polyval(fitAllB, xFitB);

wi=10;
hi=4;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

colormap jet

s1=subplot(1,2,1);

hold on
surf(edges{1},edges{2},log10(N'),'edgecolor','none')
view(2)

%axis equal
xlim([0,7])
ylim([-100,800])
caxis([0 3.5])

grid on
xlabel('Velocity (m s^{-1})');
ylabel('Offset (m)');
title(['(a) VEL above melting layer vs offset'])

plot(xFitA, yFitA,'-k','linewidth',2);
s1.SortMethod='childorder';

text(0.8,770,['Offset = ',num2str(fitAllA(1),3),' VEL - ',num2str(abs(fitAllA(2)),2)],'FontSize',12,'FontWeight','bold');

s2=subplot(1,2,2);

hold on
surf(edges{1},edges{2},log10(N2'),'edgecolor','none')
%surf(edges{1},edges{2},N2','edgecolor','none')
view(2)

xlim([0,7])
ylim([-100,800])
caxis([0 3.5])
hcb=colorbar;
hcb.Title.String='log_{10}(N)';

grid on
xlabel('Velocity (m s^{-1})');
ylabel('Offset (m)');
title(['(b) VEL below melting layer vs offset'])

plot(xFitB, yFitB,'-k','linewidth',2);
s2.SortMethod='childorder';

text(0.8,770,['Offset = ',num2str(fitAllB(1),2),' VEL - ',num2str(abs(fitAllB(2)),2)],'FontSize',12,'FontWeight','bold');

s2.Position=[0.54 0.14 0.39 0.788];
hcb.Position=[0.946 0.14 0.02 0.76];
s1.Position=[0.07 0.14 0.39 0.788];

set(gcf,'PaperPositionMode','auto')
print([figdir,'offsetVSvel'],'-dpng','-r0');
