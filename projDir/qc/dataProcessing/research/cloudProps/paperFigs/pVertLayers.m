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

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(projects)
        if ~(strcmp(projects{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(projects{jj})=in.(projects{jj}).(plotVars{ii});
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

projects=fields(plotV.cloudTopAll);

% (a) Cloud top convYoung
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
    thisVar=plotV.cloudTopAll.(projects{ii}).ConvYoungShallow;
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
    thisVar=plotV.cloudTopAll.(projects{ii}).ConvYoungMid;
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
    thisVar=plotV.cloudTopAll.(projects{ii}).ConvYoungDeep;
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

text(1:1:9,zeros(1,9)-0.75,num2str(nums),'HorizontalAlignment','center')
ylim([-1.5,15]);
xlim([0.5,9.5]);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticks(2:3:8);
xticklabels({'ConvYoungShallow','ConvYoungMid','ConvYoungDeep'});

plot([3.5,3.5],[-100,100],'-k','LineWidth',0.4);
plot([6.5,6.5],[-100,100],'-k','LineWidth',0.4);

ylabel('Cloud top (km)');
title('(a) Cloud top (km)');

text(0.7,14,'CSET','FontSize',12,'Color','g','FontWeight','bold')
text(0.7,13,'SOCRATES','FontSize',12,'Color','b','FontWeight','bold')
text(0.7,12,'OTREC','FontSize',12,'Color','r','FontWeight','bold')

% (c) Cloud top shallow CSET

load coastlines

s3=subplot(3,2,3);

hold on

thisLons=plotV.lonAll.cset.ConvYoungShallow;
thisLons=cat(1,thisLons,plotV.lonAll.cset.ConvMatureShallow);
thisLons=cat(1,thisLons,plotV.lonAll.cset.StratShallow);
thisLats=plotV.latAll.cset.ConvYoungShallow;
thisLats=cat(1,thisLats,plotV.latAll.cset.ConvMatureShallow);
thisLats=cat(1,thisLats,plotV.latAll.cset.StratShallow);

plotVar=plotV.cloudTopAll.cset.ConvYoungShallow;
plotVar=cat(1,plotVar,plotV.cloudTopAll.cset.ConvMatureShallow);
plotVar=cat(1,plotVar,plotV.cloudTopAll.cset.StratShallow);

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 4]);
s3.Colormap=jet;
cb3=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(c) CSET: cloud top shallow clouds (km)');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (e) Cloud top shallow OTREC

s5=subplot(3,2,5);

hold on

thisLons=plotV.lonAll.otrec.ConvYoungShallow;
thisLons=cat(1,thisLons,plotV.lonAll.otrec.ConvMatureShallow);
thisLons=cat(1,thisLons,plotV.lonAll.otrec.StratShallow);
thisLats=plotV.latAll.otrec.ConvYoungShallow;
thisLats=cat(1,thisLats,plotV.latAll.otrec.ConvMatureShallow);
thisLats=cat(1,thisLats,plotV.latAll.otrec.StratShallow);

plotVar=plotV.cloudTopAll.otrec.ConvYoungShallow;
plotVar=cat(1,plotVar,plotV.cloudTopAll.otrec.ConvMatureShallow);
plotVar=cat(1,plotVar,plotV.cloudTopAll.otrec.StratShallow);

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([1 5]);
s5.Colormap=jet;
cb5=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(e) OTREC: cloud top shallow clouds (km)');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on


% Layers
edges=0.5:1:9.5;

% (b) Shallow all

countsPerc=[];
medians=[];
nums=[];
for ii=1:length(projects)
    thisVar=cat(1,plotV.cloudLayersAll.(projects{ii}).ConvYoungShallow,...
        plotV.cloudLayersAll.(projects{ii}).ConvMatureShallow,...
        plotV.cloudLayersAll.(projects{ii}).StratShallow);
    thisNum=sum(~isnan(thisVar));
    counts=histcounts(thisVar,edges);
    if thisNum<10
        counts(:)=0;
    end
    countsPerc=cat(1,countsPerc,counts./thisNum.*100);
    medians=cat(1,medians,median(thisVar,'omitnan'));
    nums=cat(1,nums,thisNum);
end

s2=subplot(4,2,2);
b=bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor','flat');
for kk=1:size(countsPerc,1)
    b(kk).CData=cols(kk,:);
end

grid on
box on

xlim([edges(1),edges(end)]);
ylim([0,100]);

ylabel('Percent (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

title('(b) Cloud layers shallow clouds');

xlabel('Number of cloud layers');
legend({'CSET','SOCRATES','OTREC'});

% (d1) Mid SOCRATES

countsPerc=[];
medians=[];
nums=[];

midClouds={'ConvYoungMid','ConvMatureMid','StratMid'};
for ii=1:length(midClouds)
    thisVar=plotV.cloudLayersAll.socrates.(midClouds{ii});
    thisNum=sum(~isnan(thisVar));
    counts=histcounts(thisVar,edges);
    if thisNum<10
        counts(:)=0;
    end
    countsPerc=cat(1,countsPerc,counts./thisNum.*100);
    medians=cat(1,medians,median(thisVar,'omitnan'));
    nums=cat(1,nums,thisNum);
end

colsSoc=[0,0,0.5;0,0.5,1;0,1,1];

s41=subplot(4,2,4);
b=bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor','flat');
for kk=1:size(countsPerc,1)
    b(kk).CData=colsSoc(kk,:);
end

grid on
box on

xlim([edges(1),edges(end)]);
ylim([0,100]);

ylabel('Percent (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

title('(d) Cloud layers mid-level clouds');
legend({'SOCRATES ConvYoungMid','SOCRATES ConvMatureMid','SOCRATES StratMid'});

% (d2) Mid OTREC

countsPerc=[];
medians=[];
nums=[];

midClouds={'ConvYoungMid','ConvMatureMid','StratMid'};
for ii=1:length(midClouds)
    thisVar=plotV.cloudLayersAll.otrec.(midClouds{ii});
    thisNum=sum(~isnan(thisVar));
    counts=histcounts(thisVar,edges);
    if thisNum<10
        counts(:)=0;
    end
    countsPerc=cat(1,countsPerc,counts./thisNum.*100);
    medians=cat(1,medians,median(thisVar,'omitnan'));
    nums=cat(1,nums,thisNum);
end

colsOt=[0.5,0,0;1,0,0.2;1,0.5,1];

s42=subplot(4,2,6);
b=bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor','flat');
for kk=1:size(countsPerc,1)
    b(kk).CData=colsOt(kk,:);
end

grid on
box on

xlim([edges(1),edges(end)]);
ylim([0,60]);

ylabel('Percent (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

xlabel('Number of cloud layers');
legend({'OTREC ConvYoungMid','OTREC ConvMatureMid','OTREC StratMid'});

% (f) Deep OTREC

countsPerc=[];
medians=[];
nums=[];

deepClouds={'ConvYoungDeep','ConvMatureDeep','StratDeep'};
for ii=1:length(deepClouds)
    thisVar=plotV.cloudLayersAll.otrec.(deepClouds{ii});
    thisNum=sum(~isnan(thisVar));
    counts=histcounts(thisVar,edges);
    if thisNum<10
        counts(:)=0;
    end
    countsPerc=cat(1,countsPerc,counts./thisNum.*100);
    medians=cat(1,medians,median(thisVar,'omitnan'));
    nums=cat(1,nums,thisNum);
end

colsOt=[0.5,0,0;1,0,0.2;1,0.5,1];

s6=subplot(4,2,8);
b=bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor','flat');
for kk=1:size(countsPerc,1)
    b(kk).CData=colsOt(kk,:);
end

grid on
box on

xlim([edges(1),edges(end)]);
ylim([0,60]);

ylabel('Percent (%)')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

title('(f) Cloud layers deep clouds');
xlabel('Number of cloud layers');
legend({'OTREC ConvYoungDeep','OTREC ConvMatureDeep','OTREC StratDeep'});

s1.Position=[0.045,0.71,0.4,0.265];
s2.Position=[0.55,0.71,0.44,0.265];
s3.Position=[0.045,0.38,0.4,0.265];
s41.Position=[0.55,0.515,0.44,0.13];
s42.Position=[0.55,0.38,0.44,0.13];
s5.Position=[0.045,0.045,0.4,0.265];
s6.Position=[0.55,0.045,0.44,0.265];

cb3.Position=[0.455,0.38,0.015,0.265];
cb5.Position=[0.455,0.045,0.015,0.265];

s2.YTick=0:20:80;
s41.YTick=0:20:80;
s42.YTick=0:20:40;

set(gcf,'PaperPositionMode','auto')
print([figdir,'pVertLayers.png'],'-dpng','-r0');
