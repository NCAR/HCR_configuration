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

xTickLoc=[2,5,8,10.5,12.5,14,15,16];
xLab={'ConvShallow','ConvStratShallow','StratShallow','ConvMid','StratMid', ...
    'CD','CSD','SD'};
colsAll=[0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0; ...
    0,0,1;1,0,0;0,0,1;1,0,0;1,0,0;1,0,0;1,0,0];
projectsAll={'ConvShallowC','ConvShallowS','ConvShallowO', ...
    'ConvStratShallowC','ConvStratShallowS','ConvStratShallowO','StratShallowC','StratShallowS','StratShallowO', ...
    'ConvMidS','ConvMidO','StratMidS','StratMidO', ...
    'ConvDeep','ConvStratDeep','StratDeep'};

close all

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on','renderer','painters');

% (a) Fraction
varName='frac';
[allVars,groups,nums]=getVarShafts(plotV,varName);

s1=subplot(4,1,1);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,16.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.1,1]);
text(1:1:16,zeros(1,16)-0.05,num2str(nums),'HorizontalAlignment','center');
yticks(0:0.25:1);
ylabel('Fraction');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);

title('(a) Precipitation fraction');

% (b) Length
varName='shaftKM';
[allVars,groups,nums]=getVarShafts(plotV,varName);

s2=subplot(4,1,2);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,16.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-25,250]);
text(1:1:16,zeros(1,16)-12.5,num2str(nums),'HorizontalAlignment','center');
yticks(0:50:250);
ylabel('Width (km)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);

title('(b) Shaft width (km)');

% (c) Mean reflecitivy
varName='meanRef';
[allVars,groups,nums]=getVarShafts(plotV,varName);

s3=subplot(4,1,3);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,16.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-45,10]);
text(1:1:16,zeros(1,16)-42.5,num2str(nums),'HorizontalAlignment','center');
yticks(-40:10:10);
ylabel('Reflectivity (dBZ)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);

title('(c) Mean reflectivity (dBZ)');

% (d) Mean Vel
varName='meanVel';
[allVars,groups,nums]=getVarShafts(plotV,varName);

s4=subplot(4,1,4);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,16.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-1.25,5]);
text(1:1:16,zeros(1,16)-1,num2str(nums),'HorizontalAlignment','center');
yticks(0:10);
ylabel('Velocity (m s^{-1})');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([3.5,3.5],[-100,1000],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([14.5,14.5],[-100,1000],'-k','LineWidth',0.4);

title('(d) Mean velocity (m s^{-1})');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.2);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');

s1.Position=[0.05,0.77,0.945,0.19];
s2.Position=[0.05,0.525,0.945,0.19];
s3.Position=[0.05,0.28,0.945,0.19];
s4.Position=[0.05,0.035,0.945,0.19];

set(gcf,'PaperPositionMode','auto')
print([figdir,'shafts.png'],'-dpng','-r0');
