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

xTickLoc=[1.5,4,6.5,8.5,10.5,12.5,14.5,17,19,20,21];
xLab={'CloudLow','ConvShallow','ConvStratShallow','StratShallow','CloudMid','ConvMid','StratMid', ...
    'CloudHigh','CD','CSD','SD'};
colsAll=[0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;0,1,0;0,0,1;0,1,0;0,0,1; ...
    0,0,1;1,0,0;0,0,1;1,0,0;0,0,1;1,0,0;0,1,0;0,0,1;1,0,0;1,0,0;1,0,0;1,0,0];
projectsAll={'CloudLowS','CloudLowO','ConvShallowC','ConvShallowS','ConvShallowO', ...
    'ConvStratShallowC','ConvStratShallowS','StratShallowC','StratShallowS', ...
    'cloudMidS','cloudMidO','ConvMidS','ConvMidO','StratMidS','StratMidO', ...
    'CloudHighC','CloudHighS','CloudHighO','ConvDeep','ConvStratDeep','StratDeep'};

close all

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1000],'visible','on','renderer','painters');

% (a) Alt percentile
varName='cloudAltPerc';
[allVars,groups,nums]=getVarUpRegs(plotV,varName);

s1=subplot(5,1,1);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,21.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-15,100]);
text(1:1:21,zeros(1,21)-7.5,num2str(nums),'HorizontalAlignment','center');
yticks(0:25:100);
ylabel('Percent (%)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([2.5,2.5],[-100,1000],'-k','LineWidth',0.4);
plot([5.5,5.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([19.5,19.5],[-100,1000],'-k','LineWidth',0.4);
plot([20.5,20.5],[-100,1000],'-k','LineWidth',0.4);

title('(a) Altitude percentile');

% (b) Area
varName='area';
[allVars,groups,nums]=getVarUpRegs(plotV,varName);

s2=subplot(5,1,2);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,21.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.2,1.4]);
text(1:1:21,zeros(1,21)-0.1,num2str(nums),'HorizontalAlignment','center');
yticks(0:0.2:1.4);
ylabel('Area (km^2)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([2.5,2.5],[-100,1000],'-k','LineWidth',0.4);
plot([5.5,5.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([19.5,19.5],[-100,1000],'-k','LineWidth',0.4);
plot([20.5,20.5],[-100,1000],'-k','LineWidth',0.4);

title('(b) Area (km^2)');

% (c) Width
varName='width';
[allVars,groups,nums]=getVarUpRegs(plotV,varName);

s3=subplot(5,1,3);
hold on
boxplot(allVars,groups,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,21.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.5,4]);
text(1:1:21,zeros(1,21)-0.25,num2str(nums),'HorizontalAlignment','center');
yticks(0:4);
ylabel('Width (km)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([2.5,2.5],[-100,1000],'-k','LineWidth',0.4);
plot([5.5,5.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([19.5,19.5],[-100,1000],'-k','LineWidth',0.4);
plot([20.5,20.5],[-100,1000],'-k','LineWidth',0.4);

title('(c) Width (km)');

% (d) Depth
varName='depth';
[allVars2,groups2,nums2]=getVarUpRegs(plotV,varName);

s4=subplot(5,1,4);
hold on
boxplot(allVars2,groups2,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,21.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.3,2]);
text(1:1:21,zeros(1,21)-0.15,num2str(nums2),'HorizontalAlignment','center');
yticks(0:0.5:2);
ylabel('Depth (km)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([2.5,2.5],[-100,1000],'-k','LineWidth',0.4);
plot([5.5,5.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([19.5,19.5],[-100,1000],'-k','LineWidth',0.4);
plot([20.5,20.5],[-100,1000],'-k','LineWidth',0.4);

title('(d) Depth (km)');

% (e) Depth/width
varName='depth';
[allVars2,groups2,nums2]=getVarUpRegs(plotV,varName);

s5=subplot(5,1,5);
hold on
boxplot(allVars2./allVars,groups2,'ColorGroup',projectsAll,'Colors',colsAll,'Whisker',1.5,'Symbol','');

xlim([0.5,21.5]);
xticks(xTickLoc);
xticklabels(xLab);
xtickangle(0)
set(gca,'TickLength',[0 .01]);

ylim([-0.3,2.2]);
text(1:1:21,zeros(1,21)-0.15,num2str(nums2),'HorizontalAlignment','center');
yticks(0:0.5:2);
ylabel('Depth fraction');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

plot([2.5,2.5],[-100,1000],'-k','LineWidth',0.4);
plot([5.5,5.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([9.5,9.5],[-100,1000],'-k','LineWidth',0.4);
plot([11.5,11.5],[-100,1000],'-k','LineWidth',0.4);
plot([13.5,13.5],[-100,1000],'-k','LineWidth',0.4);
plot([15.5,15.5],[-100,1000],'-k','LineWidth',0.4);
plot([18.5,18.5],[-100,1000],'-k','LineWidth',0.4);
plot([19.5,19.5],[-100,1000],'-k','LineWidth',0.4);
plot([20.5,20.5],[-100,1000],'-k','LineWidth',0.4);

title('(e) Depth/width');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.2);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.2,'LineStyle','-','Color','k');

s1.Position=[0.045,0.825,0.945,0.15];
s2.Position=[0.045,0.625,0.945,0.15];
s3.Position=[0.045,0.425,0.945,0.15];
s4.Position=[0.045,0.225,0.945,0.15];
s5.Position=[0.045,0.025,0.945,0.15];

set(gcf,'PaperPositionMode','auto')
print([figdir,'upRegs.png'],'-dpng','-r0');
