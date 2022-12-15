% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

for ii=1:length(project)
    if strcmp(project{ii},'cset')
        qcVersion='v3.0';
        quality='qc3';
    else
        qcVersion='v3.1';
        quality='qc3';
    end

    cfDir=HCRdir(project{ii},quality,qcVersion,freqData);

    indir=[cfDir(1:end-5),'cloudProps/'];

    in.(project{ii})=load([indir,project{ii},'_cloudProps.mat']);
end

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(project)
        if ~(strcmp(project{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(project{jj})=in.(project{jj}).(plotVars{ii});
        end
    end
end

%% Plot

disp('Plotting ...')


%% Plot

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1100],'visible','on','renderer','painters');

cols=[0,1,0; ...
    0,0,1; ...
    1,0,0];

projects=fields(plotV.cloudBaseAll);

% (a) Cloud base CloudLow
allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.cloudBaseAll.(projects{ii}).CloudLow;
    thisNum=sum(~isnan(thisVar));
    if thisNum<20
        thisVar=nan;
    end
    nums=cat(1,nums,thisNum);
    allVars=cat(1,allVars,thisVar);
    groups=cat(1,groups,repmat({projects{ii}},length(thisVar),1));
    medians=cat(1,medians,median(thisVar,'omitnan'));
end

s1=subplot(3,2,1);
boxplot(allVars,groups,'ColorGroup',projects,'Colors',cols,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1.25:1:3.25,medians,num2str(nums))
ylim([0,4.5]);
yticks(0:5);

xticklabels({'CSET','SOCRATES','OTREC'});

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel('Cloud base (km)');
title('(a) Cloud base CloudLow (km)');

% (b) Cloud base OTREC

load coastlines

s2=subplot(3,2,2);

hold on

thisLons=plotV.lonAll.otrec.CloudLow;
thisLats=plotV.latAll.otrec.CloudLow;

scatter(thisLons,thisLats,25,plotV.cloudBaseAll.otrec.CloudLow,'filled');

if length(plotV.cloudBaseAll.otrec.CloudLow)>15
    perc=prctile(plotV.cloudBaseAll.otrec.CloudLow,[5,95]);
    if length(unique(perc))>1
        caxis(perc);
    end
end
caxis([0 4]);
s2.Colormap=jet;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(b) OTREC: cloud base CloudLow (km)');

ylabel('Latitude (deg)')

grid on
box on

% (c) Cloud depth all cloud
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

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.cloudDepthAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.cloudDepthAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.cloudDepthAll.(projects{ii}).CloudHigh;
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

text(0.65:1:8.65,zeros(1,9)-0.25,num2str(nums))
ylim([-0.5,5]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-2,6],'-k','LineWidth',0.4);
plot([6.5,6.5],[-2,6],'-k','LineWidth',0.4);

ylabel('Cloud depth (km)');
title('(c) Cloud depth (km)');

text(0.7,4.8,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,4.5,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,4.2,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (d) Cloud base OTREC

load coastlines

s4=subplot(3,2,4);

hold on

thisLons=plotV.lonAll.otrec.CloudLow;
thisLats=plotV.latAll.otrec.CloudLow;

scatter(thisLons,thisLats,25,plotV.cloudDepthAll.otrec.CloudLow,'filled');

if length(plotV.cloudBaseAll.otrec.CloudLow)>15
    perc=prctile(plotV.cloudBaseAll.otrec.CloudLow,[5,95]);
    if length(unique(perc))>1
        caxis(perc);
    end
end
caxis([0 4]);
s4.Colormap=jet;
cb=colorbar;

xlim(lonLims(ii,:));
ylim(latLims(ii,:));

plot(coastlon,coastlat,'-k')
title('(d) OTREC: cloud depth CloudLow (km)');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (e) Cloud length all cloud
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

allVars=[];
groups={};
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=plotV.cloudLengthAll.(projects{ii}).CloudLow;
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
    thisVar=plotV.cloudLengthAll.(projects{ii}).CloudMid;
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
    thisVar=plotV.cloudLengthAll.(projects{ii}).CloudHigh;
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

text(0.65:1:8.65,zeros(1,9)-5,num2str(nums))
ylim([-10,90]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'CloudLow','CloudMid','CloudHigh'});

plot([3.5,3.5],[-20,160],'-k','LineWidth',0.4);
plot([6.5,6.5],[-20,160],'-k','LineWidth',0.4);

ylabel('Cloud width (km)');
title('(e) Cloud width (km)');

text(0.7,86,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,81,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,76,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (f) Cloud length CSET

load coastlines

s6=subplot(3,2,6);

hold on

thisLons=plotV.lonAll.cset.CloudLow;
thisLats=plotV.latAll.cset.CloudLow;

scatter(thisLons,thisLats,25,plotV.cloudLengthAll.cset.CloudLow,'filled');

caxis([0 20]);
s6.Colormap=jet;
cb2=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(f) CSET: cloud width CloudLow (km)');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.05,0.705,0.42,0.27];
s2.Position=[0.525,0.705,0.42,0.27];
s3.Position=[0.05,0.385,0.42,0.27];
s4.Position=[0.525,0.385,0.42,0.27];
s5.Position=[0.05,0.045,0.42,0.27];
s6.Position=[0.525,0.045,0.42,0.27];

cb.Position=[0.955,0.385,0.02,0.59];
cb2.Position=[0.955,0.045,0.02,0.27];

set(gcf,'PaperPositionMode','auto')
print([figdir,'nonpDim.png'],'-dpng','-r0');
