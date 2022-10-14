function plotStats(propIn,edges,xlab,figname,classTypes,colmapCC)

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'visible','on');

for jj=1:length(classTypes)
    thisVar=propIn.(classTypes{jj});
    counts=histcounts(thisVar,edges);
    countsPerc=counts./sum(~isnan(thisVar)).*100;

    subplot(4,3,jj);
    % Get x-axis right
    edgeSpace=edges(3)-edges(2);
    plotEdges=[edges(2)-edgeSpace,edges(2:end-1),edges(end-1)+edgeSpace];
    b=bar(plotEdges(1:end-1)+edgeSpace/2,countsPerc,1,'FaceColor','flat');
    for kk=1:size(countsPerc,1)
        b(kk).CData=cols(kk,:);
    end

    grid on
    box on

    xlim([plotEdges(1),plotEdges(end)]);
    xlabel(xlab);
    ylabel('Percent (%)')

    title([classTypes{jj},' (',num2str(sum(~isnan(thisVar))),'). Median: ',num2str(median(thisVar,'omitnan'),'%.2f')]);
end

sgtitle(xlab,'FontSize',14,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print(figname,'-dpng','-r0');

end