% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

projects={'cset';'socrates';'otrec'}; %socrates, aristo, cset, otrec
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

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(projects)
        if ~(strcmp(projects{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(projects{jj})=in.(projects{jj}).(plotVars{ii});
        end
    end
end

%% Plot

cols3=[0,1,0; ...
    0,0,1; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,0,0];

projects3=cat(1,projects,projects,projects);

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1100],'visible','on','renderer','painters');

% (a) Mean convectivity ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvYoungShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvYoungMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvYoungDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s1=subplot(3,3,1);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Convectivity');
title('(a) Mean convectivity ConvYoung');

% (b) Mean convectivity ConvMature

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvMatureShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvMatureMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).ConvMatureDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s2=subplot(3,3,2);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s2.YTickLabel=[];
title('(b) Mean convectivity ConvMature');

% (c) Mean convectivity Strat

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).StratShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).StratMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).StratDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s3=subplot(3,3,3);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s3.YTickLabel=[];
title('(c) Mean convectivity Strat');

text(0.7,0.96,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,0.9,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,0.84,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (d) Max convectivity ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvYoungShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvYoungMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvYoungDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s4=subplot(3,3,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
%set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','');
%set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
%set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
%set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
%set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Convectivity');
title('(d) Maximum convectivity ConvYoung');
hold off

% (e) Max convectivity ConvMature

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvMatureShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvMatureMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).ConvMatureDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s5=subplot(3,3,5);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
% set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
% set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
% set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s5.YTickLabel=[];
title('(e) Maximum convectivity ConvMature');

% (f) Max convectivity Strat

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).StratShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).StratMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).StratDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s6=subplot(3,3,6);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
% set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
% set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
% set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center');
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s6.YTickLabel=[];
title('(f) Maximum convectivity Strat');

% (g) Max reflectivity ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvYoungShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvYoungMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvYoungDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s7=subplot(3,3,7);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
% set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
% set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
% set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-32.5,num2str(nums),'HorizontalAlignment','center');
ylim([-35,30]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Reflectivity (dBZ)');
title('(g) Maximum reflectivity ConvYoung');

% (h) Max reflectivity ConvMature

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvMatureShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvMatureMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).ConvMatureDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s8=subplot(3,3,8);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
% set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
% set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
% set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
% set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-32.5,num2str(nums),'HorizontalAlignment','center');
ylim([-35,30]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

s8.YTickLabel=[];
title('(h) Maximum reflectivity ConvMature');

% (i) Max reflectivity Strat

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).StratShallow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).StratMid;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).StratDeep;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s9=subplot(3,3,9);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx1=findobj('Tag','boxplot');
set(findobj(bx1,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx1,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx1,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx1,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx1,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx1,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-32.5,num2str(nums),'HorizontalAlignment','center');
ylim([-35,30]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

s9.YTickLabel=[];
title('(i) Maximum reflectivity Strat');

lines = s4.Children;
box4m=findobj(bx1,'Tag','Median');
box4b=findobj(bx1,'Tag','Box');
box4m(54).Color=[0,1,0];
box4b(54).Color=[0,1,0];
box4m(53).Color=[0,0,1];
box4b(53).Color=[0,0,1];
box4m(52).Color=[1,0,0];
box4b(52).Color=[1,0,0];
box4m(50).Color=[0,0,1];
box4b(50).Color=[0,0,1];
box4m(49).Color=[1,0,0];
box4b(49).Color=[1,0,0];
box4m(46).Color=[1,0,0];
box4b(46).Color=[1,0,0];

s1.Position=[0.05,0.69,0.31,0.28];
s2.Position=[0.365,0.69,0.31,0.28];
s3.Position=[0.68,0.69,0.31,0.28];

s4.Position=[0.05,0.36,0.31,0.28];
s5.Position=[0.365,0.36,0.31,0.28];
s6.Position=[0.68,0.36,0.31,0.28];

s7.Position=[0.05,0.03,0.31,0.28];
s8.Position=[0.365,0.03,0.31,0.28];
s9.Position=[0.68,0.03,0.31,0.28];

set(gcf,'PaperPositionMode','auto')
print([figdir,'pConv.png'],'-dpng','-r0');
