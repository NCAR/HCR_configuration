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


fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on');

cols=[0,1,0; ...
    0,0,1; ...
    1,0,0];

projects=fields(plotV.cloudBaseAll);

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

s=subplot(4,3,jj);
boxplot(allVars,groups,'ColorGroup',projects,'Colors',cols,'Whisker',1.5,'Symbol','+k');

bx=findobj('Tag','boxplot');
set(findobj(bx,'Tag','Box'),'LineWidth',1.5);
set(findobj(bx,'Tag','Median'),'LineWidth',1.5,'Color','k');
set(findobj(bx,'Tag','Upper Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Whisker'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Upper Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');
set(findobj(bx,'Tag','Lower Adjacent Value'),'LineWidth',1.5,'LineStyle','-','Color','k');

text(1.25:1:3.25,medians,num2str(nums))
ylim(ylimits)

set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel(ylab);
title([classTypes{jj}],'Color',colmapCC(jj,:));

set(gcf,'PaperPositionMode','auto')
print([figname,'_box.png'],'-dpng','-r0');
