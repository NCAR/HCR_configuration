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

% (a) Mean up vel ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvYoungShallow;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvYoungMid;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvYoungDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Velocity (m s^{-1})');
title('(a) Mean upward motion ConvYoung');

% (b) Mean convectivity ConvMature

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvMatureShallow;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvMatureMid;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).ConvMatureDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s2.YTickLabel=[];
title('(b) Mean upward motion ConvMature');

% (c) Mean up vel Strat

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).StratShallow;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).StratMid;
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
    thisVar=plotV.upMeanStrengthAll.(projects{ii}).StratDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-0.1,num2str(nums),'HorizontalAlignment','center');
ylim([-0.2,3]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

s3.YTickLabel=[];
title('(c) Mean upward motion Strat');

text(0.7,2.8,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,2.65,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,2.5,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (d) Max up motion ConvYoung

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvYoungShallow;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvYoungMid;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvYoungDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-1.5,num2str(nums),'HorizontalAlignment','center');
ylim([-3,25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Velocity (m s^{-1})');
title('(d) Maximum upward motion ConvYoung');

% (e) Max up motion ConvMature

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvMatureShallow;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvMatureMid;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).ConvMatureDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-1.5,num2str(nums),'HorizontalAlignment','center');
ylim([-3,25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

s5.YTickLabel=[];
title('(e) Maximum upward motion ConvMature');

% (f) Max convectivity Strat

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).StratShallow;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).StratMid;
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
    thisVar=plotV.upMaxStrengthAll.(projects{ii}).StratDeep;
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
boxplot(allVars,groups,'ColorGroup',projects3,'Colors',cols3,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1:1:9,zeros(1,9)-1.5,num2str(nums),'HorizontalAlignment','center');
ylim([-3,25]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'Shallow','Mid','Deep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

s6.YTickLabel=[];
title('(f) Maximum upward motion Strat');

% (g) Mean up vel conv young shallow OTREC

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

load coastlines

s7=subplot(3,2,5);

hold on

thisLons=plotV.lonAll.otrec.ConvYoungShallow;
thisLats=plotV.latAll.otrec.ConvYoungShallow;

plotVar=plotV.upMeanStrengthAll.otrec.ConvYoungShallow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 1.5]);
colmap=jet;
s7.Colormap=colmap(10:230,:);
cb7=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(g) OTREC: mean upward motion ConvYoungShallow (m s^{-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (h) Mean up vel conv young mid OTREC

s8=subplot(3,2,6);

hold on

thisLons=plotV.lonAll.otrec.ConvYoungMid;
thisLats=plotV.latAll.otrec.ConvYoungMid;

plotVar=plotV.upMeanStrengthAll.otrec.ConvYoungMid;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 2.5]);
colmap=jet;
s8.Colormap=colmap(10:230,:);
cb8=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(h) OTREC: mean upward motion ConvYoungMid (m s^{-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.05,0.7,0.31,0.28];
s2.Position=[0.365,0.7,0.31,0.28];
s3.Position=[0.68,0.7,0.31,0.28];

s4.Position=[0.05,0.375,0.31,0.28];
s5.Position=[0.365,0.375,0.31,0.28];
s6.Position=[0.68,0.375,0.31,0.28];

s7.Position=[0.05,0.045,0.38,0.28];
s8.Position=[0.55,0.045,0.38,0.28];

set(gcf,'PaperPositionMode','auto')
print([figdir,'pUpVel.png'],'-dpng','-r0');
