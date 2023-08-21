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


%% Plot

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on','renderer','painters');

load coastlines

% (a) Cloud base CloudLow

s1=subplot(2,2,1);

hold on

thisLons=plotV.lonAll.otrec.CloudLow;

thisLats=plotV.latAll.otrec.CloudLow;

plotVar=plotV.cloudBaseAll.otrec.CloudLow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 4]);
s1.Colormap=jet;
cb1=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(a) Cloud base (km). Low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (b) Cloud depth CloudLow

s2=subplot(2,2,2);

hold on

thisLons=plotV.lonAll.otrec.CloudLow;

thisLats=plotV.latAll.otrec.CloudLow;

plotVar=plotV.cloudDepthAll.otrec.CloudLow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 5]);
s2.Colormap=jet;
cb2=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(b) Cloud depth (km). Low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (c) Cloud top

s3=subplot(2,2,3);

hold on

thisLons=cat(1,plotV.lonAll.otrec.CloudLow,plotV.lonAll.otrec.ConvYoungShallow, ...
    plotV.lonAll.otrec.ConvMatureShallow,plotV.lonAll.otrec.StratShallow);

thisLats=cat(1,plotV.latAll.otrec.CloudLow,plotV.latAll.otrec.ConvYoungShallow, ...
    plotV.latAll.otrec.ConvMatureShallow,plotV.latAll.otrec.StratShallow);

plotVar=cat(1,plotV.cloudTopAll.otrec.CloudLow,plotV.cloudTopAll.otrec.ConvYoungShallow, ...
    plotV.cloudTopAll.otrec.ConvMatureShallow,plotV.cloudTopAll.otrec.StratShallow);

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([1 5]);
s3.Colormap=jet;
cb3=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(c) Cloud top (km). Shallow and low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (d) Width all

allVars=[];
groups={};
medians=[];
nums=[];


thisVar=plotV.cloudLengthAll.otrec.CloudLow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudLow'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.ConvYoungShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvShallow'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.ConvMatureShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratShallow'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.StratShallow;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratShallow'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.CloudMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudMid'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.ConvYoungMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvMid'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.StratMid;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratMid'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

s4=subplot(2,2,4);
hold on
boxplot(allVars,groups,'Whisker',1.5,'Symbol','');
%boxplot(allVars,groups,'ColorGroup',projects2,'Colors',cols2,'Whisker',1.5,'Symbol','');

text(1:1:7,zeros(1,7)-5,num2str(nums),'HorizontalAlignment','center');
ylim([-10,100]);
xlim([0.5,11.5]);
yticks(0:30:90);

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel('Width Low/Shallow/Mid (km)');

allVars=[];
groups={};
medians=[];
nums=[];

thisVar=plotV.cloudLengthAll.otrec.CloudHigh;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'CloudHigh'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.ConvYoungDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvDeep'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.ConvMatureDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'ConvStratDeep'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));

thisVar=plotV.cloudLengthAll.otrec.StratDeep;
thisNum=sum(~isnan(thisVar));
nums=cat(1,nums,thisNum);
allVars=cat(1,allVars,thisVar);
groups=cat(1,groups,repmat({'StratDeep'},length(thisVar),1));
medians=cat(1,medians,median(thisVar,'omitnan'));


allVars=cat(1,nan,nan,nan,nan,nan,nan,nan,allVars);
groups=cat(1,{'nan1'},{'nan2'},{'nan3'},{'nan4'},{'nan5'},{'nan6'},{'nan7'},groups);

yyaxis right
boxplot(allVars,groups,'Whisker',1.5,'Symbol','');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

numstext={};
for ii=1:4
    numstext{end+1}=num2str(nums(ii));
end
   
text(1:1:11,zeros(1,11)-33,[{''};{''};{''};{''};{''};{''};{''};numstext'],'HorizontalAlignment','center');
ylim([-67,667]);
xlim([0.5,11.5]);
yticks(0:200:600);
ylabel('Width High/Deep (km)');

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

xticklabels({'CloudLow','ConvShallow','ConvStratShallow','StratShallow','cloudMid','ConvMid','StratMid', ...
    'CloudHigh','ConvDeep','ConvStratDeep','StratDeep'});

plot([4.5,4.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.5,7.5],[-100,1000],'-k','LineWidth',0.4);
plot([7.6,7.6],[-100,1000],'-k','LineWidth',0.4);

title('(d) Width');

ax = gca;
ax.YColor = 'k';

s1.Position=[0.045,0.56,0.4,0.4];
s2.Position=[0.545,0.56,0.4,0.4];
s3.Position=[0.045,0.065,0.4,0.4];
s4.Position=[0.545,0.11,0.4,0.355];

cb1.Position=[0.455,0.56,0.02,0.4];
cb2.Position=[0.955,0.56,0.02,0.4];
cb3.Position=[0.455,0.065,0.02,0.4];

set(gcf,'PaperPositionMode','auto')
print([figdir,'otrecDims.png'],'-dpng','-r0');
