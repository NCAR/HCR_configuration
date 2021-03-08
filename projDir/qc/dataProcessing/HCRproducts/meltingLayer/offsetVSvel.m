% Analyze HCR clouds

clear all;
close all;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
indir=['/home/romatsch/plots/HCR/meltingLayer/offsetVel/'];


%% Plot

close all
edges={-10:0.5:10 -100:20:1000};

N=hist3(cat(2,velAbove,offset),'Edges',edges);
%N(N==0)=nan;

N2=hist3(cat(2,velBelow,offset),'Edges',edges);
%N2(N2==0)=nan;

f1 = figure('Position',[200 500 1300 600],'DefaultAxesFontSize',12);
colormap jet

s1=subplot(1,2,1);

hold on
%surf(edges{1},edges{2},log10(N'),'edgecolor','none')
surf(edges{1},edges{2},N','edgecolor','none')
view(2)

%axis equal
xlim([0,8])
ylim([-100,700])
caxis([0 3000])
%xticks(-40:20:60);
%yticks(-40:20:60);

grid on
xlabel('Velocity above melting layer (m s^{-1})');
ylabel('Offset (m)');
title(['Velocity vs offset'])
s1pos=s1.Position;

% Regression
% fitOrth=gmregress(compTablePlot.DBZsur,compTablePlot.DBZrhi,1);
% fitAll=[fitOrth(2) fitOrth(1)];
% xFit = -40:0.1:60;
% yFit = polyval(fitAll, xFit);
% 
% plot(xFit, yFit,'-b','linewidth',2);
% 
% ax2.SortMethod='childorder';

s2=subplot(1,2,2);

hold on
%surf(edges{1},edges{2},log10(N'),'edgecolor','none')
surf(edges{1},edges{2},N2','edgecolor','none')
view(2)

xlim([0,8])
ylim([-100,700])
caxis([0 3000])
colorbar

grid on
xlabel('Velocity above melting layer (m s^{-1})');
ylabel('Offset (m)');
title(['Velocity vs offset'])
s2pos=s2.Position;
s2.Position=[s2pos(1) s1pos(2) s1pos(3) s1pos(4)];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'offsetVSvel_',project],'-dpng','-r0');

save([figdir,project,'_velOffset.mat'],'velAbove','velBelow','offset');