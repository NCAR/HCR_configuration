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

% (a) Mean convectivity

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.meanConvAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.meanConvAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.meanConvAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s1=subplot(3,2,1);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Convectivity');
title('(a) Mean convectivity');

text(0.7,0.96,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,0.9,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,0.84,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (b) Max convectivity

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxConvAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.maxConvAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.maxConvAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s2=subplot(3,2,2);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Convectivity');
title('(b) Maximum convectivity');

% (c) Mean reflectivity

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.meanReflAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.meanReflAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.meanReflAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s3=subplot(3,2,3);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-54,num2str(nums),'HorizontalAlignment','center');
ylim([-58,25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Reflectivity (dBZ)');
title('(c) Mean reflectivity');

% (d) Max reflectivity

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.maxReflAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.maxReflAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.maxReflAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s4=subplot(3,2,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-54,num2str(nums),'HorizontalAlignment','center');
ylim([-58,25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Reflectivity (dBZ)');
title('(d) Maximum reflectivity');

% (e) Up frac

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upFracAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.upFracAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.upFracAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s5=subplot(3,2,5);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.05,num2str(nums),'HorizontalAlignment','center')
ylim([-0.1,1]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Fraction');
title('(e) Area fraction of upward motion');

% (f) Max up vel

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).CloudHigh;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s6=subplot(3,2,6);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.25,num2str(nums),'HorizontalAlignment','center');
ylim([-0.5,5]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Velocity (m s^{-1})');
title('(f) Maximum upward velocity');

s1.Position=[0.05,0.695,0.44,0.28];
s2.Position=[0.55,0.695,0.44,0.28];
s3.Position=[0.05,0.365,0.44,0.28];
s4.Position=[0.55,0.365,0.44,0.28];
s5.Position=[0.05,0.03,0.44,0.28];
s6.Position=[0.55,0.03,0.44,0.28];

set(gcf,'PaperPositionMode','auto')
print([figdir,'nonpConv.png'],'-dpng','-r0');
