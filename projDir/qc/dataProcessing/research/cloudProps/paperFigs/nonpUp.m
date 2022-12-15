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


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on','renderer','painters');

% (a) Altitude percentile

allVars=[];
groups={};
medians=[];
nums=[];

varIn='cloudAltPerc';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).CloudLow;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudMid;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudHigh;
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

s1=subplot(2,2,1);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Percent (%)');
title('(a) Altitude percentile');

text(0.7,16,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,10,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,4,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (b) Area

allVars=[];
groups={};
medians=[];
nums=[];

varIn='area';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).CloudLow;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudMid;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudHigh;
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

s2=subplot(2,2,2);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,1.5]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Area (km^2)');
title('(b) Area');

% (c) Width

allVars=[];
groups={};
medians=[];
nums=[];

varIn='width';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).CloudLow;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudMid;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudHigh;
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

s3=subplot(2,2,3);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.2,num2str(nums),'HorizontalAlignment','center');
ylim([-0.4,4]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Width (km)');
title('(c) Width');

% (d) Depth

allVars=[];
groups={};
medians=[];
nums=[];

varIn='depth';
for ii=1:length(projects)
     thisTable=plotV.upRegsAll.(projects{ii}).CloudLow;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudMid;
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
     thisTable=plotV.upRegsAll.(projects{ii}).CloudHigh;
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

s4=subplot(2,2,4);
hold on
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

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
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Depth (km)');
title('(d) Depth');

s1.Position=[0.05,0.535,0.44,0.43];
s2.Position=[0.55,0.535,0.44,0.43];
s3.Position=[0.05,0.035,0.44,0.43];
s4.Position=[0.55,0.035,0.44,0.43];

set(gcf,'PaperPositionMode','auto')
print([figdir,'nonpUp.png'],'-dpng','-r0');
