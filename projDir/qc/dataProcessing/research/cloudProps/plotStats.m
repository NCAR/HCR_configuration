function plotStats(propIn,edges,xlab,figname,classTypes,colmapCC)

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'visible','on');

for jj=1:length(classTypes)
    thisVar=propIn.(classTypes{jj});
    counts=histcounts(thisVar,edges);
    countsPerc=counts./sum(~isnan(thisVar)).*100;

    subplot(4,3,jj);
    bar(edges(1:end-1)+(edges(2)-edges(1))/2,countsPerc,1,'FaceColor',colmapCC(jj,:));
    grid on
    box on

    xlim([edges(1),edges(end)]);
    xlabel(xlab);
    ylabel('Percent (%)')

    title([classTypes{jj},' (',num2str(sum(~isnan(thisVar))),'). Median: ',num2str(median(thisVar,'omitnan'),'%.2f')]);
end

sgtitle(xlab,'FontSize',14,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print(figname,'-dpng','-r0');

end