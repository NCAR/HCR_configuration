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

minlength=2.5;

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

cols2=[0,1,0; ...
    0,0,1; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,0,0];

cols=[0,1,0; ...
    0,0,1; ...
    1,0,0];

projects2=cat(1,projects,projects);
projects3=cat(1,projects,projects,projects);

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1100],'visible','on','renderer','painters');

% (a) Precipitation fraction

allVars=[];
groups={};
medians=[];
nums=[];

varIn='frac';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s1=subplot(4,3,1);
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

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Fraction');
title('(a) Precipitation fraction Conv');

text(0.7,0.17,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,0.1,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,0.03,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (b) Precipitation fraction ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s2=subplot(4,3,2);
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

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(b) Precipitation fraction ConvStrat');

% (c) Precipitation fraction strat

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s3=subplot(4,3,3);
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

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(c) Precipitation fraction Strat');

% (d) Precipitation length ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='shaftKM';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
                thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s4=subplot(4,3,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','');

text(1:1:6,zeros(1,6)-2.5,num2str(nums),'HorizontalAlignment','center');
ylim([-5,50]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
ylabel('Length (km)');

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

allVars=cat(1,nan,nan,nan,nan,nan,nan,allVars);
groups=cat(1,{'nan1'},{'nan2'},{'nan3'},{'nan4'},{'nan5'},{'nan6'},groups);

yyaxis right
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

numstext={};
for ii=1:3
    numstext{end+1}=num2str(nums(ii));
end
   
text(1:1:9,zeros(1,9)-12.5,[{''};{''};{''};{''};{''};{''};numstext'],'HorizontalAlignment','center');
ylim([-25,250]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,500],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,500],'-k','LineWidth',0.4);
plot([6.6,6.6],[-100,500],'-k','LineWidth',0.4);

title('(d) Shaft length Conv');

ax = gca;
ax.YColor = 'k';

% (e) Length ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

varIn='shaftKM';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
                thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s5=subplot(4,3,5);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','');

text(1:1:6,zeros(1,6)-2.5,num2str(nums),'HorizontalAlignment','center');
ylim([-5,50]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
ylabel('Length (km)');

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

allVars=cat(1,nan,nan,nan,nan,nan,nan,allVars);
groups=cat(1,{'nan1'},{'nan2'},{'nan3'},{'nan4'},{'nan5'},{'nan6'},groups);

yyaxis right
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

numstext={};
for ii=1:3
    numstext{end+1}=num2str(nums(ii));
end
   
text(1:1:9,zeros(1,9)-12.5,[{''};{''};{''};{''};{''};{''};numstext'],'HorizontalAlignment','center');
ylim([-25,250]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,500],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,500],'-k','LineWidth',0.4);
plot([6.6,6.6],[-100,500],'-k','LineWidth',0.4);

title('(e) Shaft length ConvStrat');

ax = gca;
ax.YColor = 'k';

% (g) Precipitation length strat

allVars=[];
groups={};
medians=[];
nums=[];

varIn='shaftKM';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
                thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'2']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s6=subplot(4,3,6);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','');

text(1:1:6,zeros(1,6)-2.5,num2str(nums),'HorizontalAlignment','center');
ylim([-5,50]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
ylabel('Length (km)');

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

allVars=cat(1,nan,nan,nan,nan,nan,nan,allVars);
groups=cat(1,{'nan1'},{'nan2'},{'nan3'},{'nan4'},{'nan5'},{'nan6'},groups);

yyaxis right
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

numstext={};
for ii=1:3
    numstext{end+1}=num2str(nums(ii));
end
   
text(1:1:9,zeros(1,9)-12.5,[{''};{''};{''};{''};{''};{''};numstext'],'HorizontalAlignment','center');
ylim([-25,250]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,500],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,500],'-k','LineWidth',0.4);
plot([6.6,6.6],[-100,500],'-k','LineWidth',0.4);

title('(f) Shaft length Strat');

ax = gca;
ax.YColor = 'k';

% (g) Mean refl convYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='meanRef';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s7=subplot(4,3,7);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-42.5,num2str(nums),'HorizontalAlignment','center');
ylim([-45,10]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Reflectivity (dBZ)');
title('(g) Mean reflectivity Conv');

% (h) Mean ref ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s8=subplot(4,3,8);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-42.5,num2str(nums),'HorizontalAlignment','center');
ylim([-45,10]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(h) Mean reflectivity ConvStrat');

% (i) Mean refl strat

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s9=subplot(4,3,9);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-42.5,num2str(nums),'HorizontalAlignment','center');
ylim([-45,10]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(i) Mean reflectivity Strat');

% (j) Mean vel convYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='meanVel';
for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s10=subplot(4,3,10);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.6,num2str(nums),'HorizontalAlignment','center');
ylim([-0.9,5]);
yticks(0:10);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Velocity (m s^{-1})');
title('(j) Mean velocity Conv');

% (k) Mean vel ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s11=subplot(4,3,11);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.6,num2str(nums),'HorizontalAlignment','center');
ylim([-0.9,5]);
yticks(0:10);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(k) Mean velocity ConvStrat');

% (l) Mean vel strat

allVars=[];
groups={};
medians=[];
nums=[];

for ii=1:length(projects)
     thisTable=plotV.precShaftsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratMid;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
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
     thisTable=plotV.precShaftsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.shaftKM<minlength,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s12=subplot(4,3,12);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.6,num2str(nums),'HorizontalAlignment','center');
ylim([-0.9,5]);
yticks(0:10);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

yticklabels('');
title('(l) Mean velocity Strat');

s1.Position=[0.055,0.774,0.31,0.205];
s2.Position=[0.37,0.774,0.31,0.205];
s3.Position=[0.685,0.774,0.31,0.205];

s4.Position=[0.1,0.524,0.235,0.205];
s5.Position=[0.415,0.524,0.235,0.205];
s6.Position=[0.73,0.524,0.235,0.205];

s7.Position=[0.055,0.274,0.31,0.205];
s8.Position=[0.37,0.274,0.31,0.205];
s9.Position=[0.685,0.274,0.31,0.205];

s10.Position=[0.055,0.024,0.31,0.205];
s11.Position=[0.37,0.024,0.31,0.205];
s12.Position=[0.685,0.024,0.31,0.205];

set(gcf,'PaperPositionMode','auto')
print([figdir,'shafts.png'],'-dpng','-r0');
