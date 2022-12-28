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

minArea=0.1;

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

% (a) Altitude percentile ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='cloudAltPerc';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-5,num2str(nums),'HorizontalAlignment','center');
ylim([-10,100]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Percent (%)');
title('(a) Altitude percentile ConvYoung');

text(0.7,20,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,12,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,4,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (b) Altitude percentile ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

varIn='cloudAltPerc';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-5,num2str(nums),'HorizontalAlignment','center');
ylim([-10,100]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(b) Altitude percentile ConvMature');

% (c) Altitude percentile Strat

allVars=[];
groups={};
medians=[];
nums=[];

varIn='cloudAltPerc';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-5,num2str(nums),'HorizontalAlignment','center');
ylim([-10,100]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(c) Altitude percentile Strat');

% (d) Area

allVars=[];
groups={};
medians=[];
nums=[];

varIn='area';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s4=subplot(4,3,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,1.25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Area (km^2)');
title('(d) Area ConvYoung');

% (e) Area ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

varIn='area';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s5=subplot(4,3,5);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,1.25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(e) Area ConvMature');

% (f) Area Strat

allVars=[];
groups={};
medians=[];
nums=[];

varIn='area';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
        thisNum=sum(~isnan(thisVar));
        if thisNum<20
            thisVar=nan;
        end
        nums=cat(1,nums,thisNum);
        allVars=cat(1,allVars,thisVar);
        groups=cat(1,groups,repmat({[projects{ii},'3']},length(thisVar),1));
        medians=cat(1,medians,median(thisVar,'omitnan'));
end

s6=subplot(4,3,6);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,1.25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(f) Area Strat');

% (g) Width ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='width';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.2,num2str(nums),'HorizontalAlignment','center');
ylim([-0.4,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Width (km)');
title('(g) Width ConvYoung');

% (h) Width ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

varIn='width';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.2,num2str(nums),'HorizontalAlignment','center');
ylim([-0.4,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(h) Width ConvMature');

% (i) Width Strat

allVars=[];
groups={};
medians=[];
nums=[];

varIn='width';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.2,num2str(nums),'HorizontalAlignment','center');
ylim([-0.4,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(i) Width Strat');

% (j) Depth ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];

varIn='depth';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvYoungDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,2]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Depth (km)');
title('(j) Depth ConvYoung');

% (k) Depth ConvMature

allVars=[];
groups={};
medians=[];
nums=[];

varIn='depth';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).ConvMatureDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,2]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(k) Depth ConvMature');

% (l) Depth Strat

allVars=[];
groups={};
medians=[];
nums=[];

varIn='depth';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).StratShallow;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratMid;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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
     thisTable=plotV.upRegsAll.(projects{ii}).StratDeep;
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        if strcmp(varIn,'cloudAltPerc')
            thisVar=abs(thisVar-100);
        end
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

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,2]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});
yticklabels('');

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

title('(l) Depth Strat');

s1.Position=[0.055,0.774,0.31,0.205];
s2.Position=[0.37,0.774,0.31,0.205];
s3.Position=[0.685,0.774,0.31,0.205];

s4.Position=[0.055,0.524,0.31,0.205];
s5.Position=[0.37,0.524,0.31,0.205];
s6.Position=[0.685,0.524,0.31,0.205];

s7.Position=[0.055,0.274,0.31,0.205];
s8.Position=[0.37,0.274,0.31,0.205];
s9.Position=[0.685,0.274,0.31,0.205];

s10.Position=[0.055,0.024,0.31,0.205];
s11.Position=[0.37,0.024,0.31,0.205];
s12.Position=[0.685,0.024,0.31,0.205];

set(gcf,'PaperPositionMode','auto')
print([figdir,'pUp.png'],'-dpng','-r0');
