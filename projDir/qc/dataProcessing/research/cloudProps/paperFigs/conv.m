% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

projects={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

for ii=1:length(projects)
    if strcmp(projects{ii},'cset')
        qcVersion='v3.0';
        quality='qc3';
    else
        qcVersion='v3.1';
        quality='qc3';
    end

    cfDir=HCRdir(projects{ii},quality,qcVersion,freqData);

    indir=[cfDir(1:end-5),'cloudProps/'];

    in.(projects{ii})=load([indir,projects{ii},'_cloudProps.mat']);
end

plotVars=fields(in.otrec);

for ii=1:length(plotVars)
    for jj=1:length(projects)
        if ~(strcmp(projects{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(projects{jj})=in.(projects{jj}).(plotVars{ii});
        end
    end
end

%% Plot

disp('Plotting ...')

xTickLoc=[2,5,8,11,13.5,15.5,17.5,20,22,23,24];
xLab={'CloudLow','ConvShallow','ConvStratShallow','StratShallow','CloudMid','ConvMid','StratMid', ...
    'CloudHigh','CD','CSD','SD'};
colsAll=[0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0; ...
    0,0,1;1,0,0;0,0,1;1,0,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;1,0,0;1,0,0;1,0,0];
projectsAll={'CloudLowC','CloudLowS','CloudLowO','ConvShallowC','ConvShallowS','ConvShallowO', ...
    'ConvStratShallowC','ConvStratShallowS','ConvStratShallowO','StratShallowC','StratShallowS','StratShallowO', ...
    'cloudMidS','cloudMidO','ConvMidS','ConvMidO','StratMidS','StratMidO', ...
    'CloudHighC','CloudHighS','CloudHighO','ConvDeep','ConvStratDeep','StratDeep'};

close all

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1100],'visible','on','renderer','painters');

% (a) Mean convectivity
varName='meanConvAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s1=subplot(7,1,1);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.2,1]);
text(1:1:24,zeros(1,24)-0.1,num2str(nums),'HorizontalAlignment','center');
yticks(0:0.25:1);
ylabel('Convectivity');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(a) Mean convectivity');

% (b) Max convectivity
varName='maxConvAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s2=subplot(7,1,2);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.2,1]);
text(1:1:24,zeros(1,24)-0.1,num2str(nums),'HorizontalAlignment','center');
yticks(0:0.25:1);
ylabel('Convectivity');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(b) Maximum convectivity');

% (c) Mean reflectivity
varName='meanReflAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s3=subplot(7,1,3);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-65,20]);
text(1:1:24,zeros(1,24)-55,num2str(nums),'HorizontalAlignment','center');
yticks(-40:20:40);
ylabel('Reflectivity (dBZ)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(c) Mean reflectivity (dBZ)');

% (d) Max reflectivity
varName='maxReflAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s4=subplot(7,1,4);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-70,40]);
text(1:1:24,zeros(1,24)-60,num2str(nums),'HorizontalAlignment','center');
yticks(-40:20:40);
ylabel('Reflectivity (dBZ)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(d) Maximum reflectivity (dBZ)');

% (e) Mean up vel
varName='upMeanStrengthAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s5=subplot(7,1,5);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.5,3]);
text(1:1:24,zeros(1,24)-0.25,num2str(nums),'HorizontalAlignment','center');
yticks(0:1:3);
ylabel('Velocity (m s^{-1})');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(e) Mean upward velocity (m s^{-1})');

% (f) Max up vel
varName='upMaxStrengthAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s6=subplot(7,1,6);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-5,24]);
text(1:1:24,zeros(1,24)-2.5,num2str(nums),'HorizontalAlignment','center');
yticks(0:6:24);
ylabel('Velocity (m s^{-1})');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(f) Maximum upward velocity (m s^{-1})');

% (g) Up frac
varName='upFracAll';
[allVars,groups,nums]=getVarAll(plotV,varName);

s7=subplot(7,1,7);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,24.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.2,0.8]);
text(1:1:24,zeros(1,24)-0.1,num2str(nums),'HorizontalAlignment','center');
yticks(0:0.2:1);
ylabel('Fraction');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([12.5,12.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);
plot([16.5,16.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([21.5,21.5],[-100,1000],'-k','LineWidth',0.4);
plot([22.5,22.5],[-100,1000],'-k','LineWidth',0.4);
plot([23.5,23.5],[-100,1000],'-k','LineWidth',0.4);

title('(g) Area fraction of upward motion');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.2);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');

s1.Position=[0.045,0.878,0.945,0.1];
s2.Position=[0.045,0.735,0.945,0.1];
s3.Position=[0.045,0.592,0.945,0.1];
s4.Position=[0.045,0.449,0.945,0.1];
s5.Position=[0.045,0.306,0.945,0.1];
s6.Position=[0.045,0.163,0.945,0.1];
s7.Position=[0.045,0.02,0.945,0.1];

set(gcf,'PaperPositionMode','auto')
print([figdir,'conv.png'],'-dpng','-r0');
