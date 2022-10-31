function plotStatsRegsProjects_box(propIn,varIn,minArea,ylab,figname,classTypes,colmapCC,ylimits)

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'visible','on');

cols=[0,1,0; ...
    0,0,1; ...
    1,0,0];

projects=fields(propIn);

for jj=1:length(classTypes)
    allVars=[];
    groups={};
    medians=[];
    nums=[];
    for ii=1:length(projects)
        thisTable=propIn.(projects{ii}).(classTypes{jj});
        thisTable(thisTable.area<minArea,:)=[];
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
end

sgtitle(ylab,'FontSize',14,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print([figname,'_box.png'],'-dpng','-r0');

end