function plotStatsRegsProjects(propIn,varIn,minArea,edges,xlab,figname,classTypes,colmapCC)

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'visible','on');

cols=[0,1,0; ...
    0,0,1; ...
    1,0.9,0];

projects=fields(propIn);

for jj=1:length(classTypes)
    countsPerc=[];
    medians=[];
    nums=[];
    for ii=1:length(projects)
        thisTable=propIn.(projects{ii}).(classTypes{jj});
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        thisNum=sum(~isnan(thisVar));
        counts=histcounts(thisVar,edges);
        if thisNum<10
            counts(:)=0;
        end
        countsPerc=cat(1,countsPerc,counts./thisNum.*100);
        medians=cat(1,medians,median(thisVar,'omitnan'));
        nums=cat(1,nums,thisNum);
    end

    s=subplot(4,3,jj);
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

    title([classTypes{jj}],'Color',colmapCC(jj,:));

    textBox=TextLocation([num2str(medians,'%.2f'),repmat(' (',length(nums),1),num2str(nums,'%.0f'),repmat(')',length(nums),1)],'Location','best');
    textPos=textBox.Position;
    textBox.Position=[textPos(1),textPos(2)-textPos(2)*0.02,textPos(3),textPos(4)];

    spos=s.Position;

    if jj==11
        legend(projects,'Location','southoutside','Orientation','horizontal');
        s.Position=spos;
    end
end

sgtitle(xlab,'FontSize',14,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print([figname,'.png'],'-dpng','-r0');

end