function plotStatsProjects(propIn,edges,xlab,figname,classTypes,colmapCC)

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');

projects=fields(propIn);

for jj=1:length(classTypes)
    countsPerc=[];
    for ii=1:length(projects)
        thisVar=propIn.(projects{ii}).(classTypes{jj});
        counts=histcounts(thisVar,edges);
        countsPerc=cat(1,countsPerc,counts./sum(~isnan(thisVar)).*100);
    end

    s=subplot(4,3,jj);
    if jj>=7
    bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1);
    else
        plot(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,'-','LineWidth',1);
    end
    grid on
    box on

    xlim([edges(1),edges(end)]);
    xlabel(xlab);
    ylabel('Percent (%)')

    title([classTypes{jj}],'Color',colmapCC(jj,:));

    spos=s.Position;

    if jj==11
        legend(projects,'Location','southoutside','Orientation','horizontal');
        s.Position=spos;
    end
end

sgtitle(xlab,'FontSize',14,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print(figname,'-dpng','-r0');

end