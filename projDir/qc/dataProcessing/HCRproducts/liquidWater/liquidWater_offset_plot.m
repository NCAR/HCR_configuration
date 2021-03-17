% Analyze HCR clouds

clear all;
close all;

project='cset';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/liquidWater/',project,'/offsets/'];

%% Load data
load([figdir,'sig0data.mat']);

%% Plot altitudes

close all
edges={-30:0.1:30 -30:0.1:30};

N=hist3(cat(2,saveAll(:,4),saveAll(:,1)),'Edges',edges);

N2=hist3(cat(2,saveAll(:,4),saveAll(:,1)-saveAll(:,3)),'Edges',edges);

% Regression
xFit = -30:0.1:30;

sig0alt=cat(2,saveAll(:,1)-saveAll(:,3),saveAll(:,4));
sig0alt(any(isnan(sig0alt),2),:)=[];

% Quadratic fit all data
fitAA=polyfit(sig0alt(:,2),sig0alt(:,1),2);
yFitAA = polyval(fitAA,xFit);

% Zoom in
sig0alt(sig0alt(:,2)>3,:)=[];
%sig0alt(sig0alt(:,1)<-5,:)=[];

N3=hist3(cat(2,sig0alt(:,2),sig0alt(:,1)),'Edges',edges);

% Linear fit
fitOrth=gmregress(sig0alt(:,2),sig0alt(:,1),1);
fitAll=[fitOrth(2) fitOrth(1)];
yFit = polyval(fitAll, xFit);

% Quadratic fit
fitA=polyfit(sig0alt(:,2),sig0alt(:,1),2);
yFitA = polyval(fitA,xFit);

wi=10;
hi=9;

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

s1=subplot(2,2,1);

hold on
surf(edges{1},edges{2},log10(N'),'edgecolor','none')
view(2)

ylim([0,20])
xlim([0,15])
caxis([0 3])

grid on
ylabel('Sig0 (db)');
xlabel('Altitude (km)');
title(['(a) Sig0 vs altitude'])

s2=subplot(2,2,2);
hold on
surf(edges{1},edges{2},log10(N'),'edgecolor','none')
view(2)

%axis equal
ylim([0,20])
xlim([0,5])
caxis([0 3])

grid on
ylabel('Sig0 (db)');
xlabel('Altitude (km)');
title(['(b) Sig0 vs altitude zoomed'])

s3=subplot(2,2,3);

hold on
surf(edges{1},edges{2},log10(N2'),'edgecolor','none')
view(2)

%axis equal
ylim([-10,10])
xlim([0,15])
caxis([0 3])

grid on
ylabel('Sig0 (db)');
xlabel('Altitude (km)');
title(['(c) Sig0meas - sig0mod vs altitude'])

l2=plot(xFit, yFitAA,'-b','linewidth',2);

s3.SortMethod='childorder';

legend([l2],{['y=',num2str(fitAA(1)),'x^{2}+',num2str(fitAA(2)),'x+',num2str(fitAA(3))]});


s4=subplot(2,2,4);
hold on
surf(edges{1},edges{2},log10(N3'),'edgecolor','none')
view(2)

%axis equal
ylim([-10,6])
xlim([0,3])
caxis([0 3])

grid on
ylabel('Sig0 (db)');
xlabel('Altitude (km)');
title(['(d) Sig0meas - sig0mod vs altitude zoomed'])

hcb=colorbar;
hcb.Title.String='log_{10}(N)';

l1=plot(xFit, yFit,'-k','linewidth',2);
l2=plot(xFit, yFitA,'-b','linewidth',2);

s4.SortMethod='childorder';

legend([l1 l2],{['y=',num2str(fitAll(1)),'x+',num2str(fitAll(2))],...
    ['y=',num2str(fitA(1)),'x^{2}+',num2str(fitA(2)),'x+',num2str(fitA(3))]});

% text(0.8,770,['Offset = ',num2str(fitAllB(1),2),' VEL - ',num2str(abs(fitAllB(2)),2)],'FontSize',12,'FontWeight','bold');

s1.Position=[0.1300    0.5838    0.3347    0.3412];
s2.Position=[0.5703    0.5838    0.3347    0.3412];
s3.Position=[0.1300    0.1100    0.3347    0.3412];
s4.Position=[0.5703    0.1100    0.3347    0.3412];
hcb.Position=[0.93   0.5838    0.0222    0.3413];

set(gcf,'PaperPositionMode','auto')
print([figdir,'sig0vsAlt'],'-dpng','-r0');
