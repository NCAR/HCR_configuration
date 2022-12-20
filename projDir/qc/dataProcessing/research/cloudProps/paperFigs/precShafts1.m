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


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on','renderer','painters');

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

s1=subplot(2,3,1);
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
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Percent (%)');
title('(a) Precipitation fraction ConvYoung');

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

s2=subplot(2,3,2);
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
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Percent (%)');
title('(b) Precipitation fraction ConvMature');

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

s3=subplot(2,3,3);
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
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Percent (%)');
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

s4=subplot(2,3,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','+k');

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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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

title('(d) Shaft length ConvYoung');

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

s5=subplot(2,3,5);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','+k');

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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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

title('(e) Shaft length ConvMature');

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

s6=subplot(2,3,6);
hold on
boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','+k');

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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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

set(gcf,'PaperPositionMode','auto')
print([figdir,'shafts1.png'],'-dpng','-r0');
