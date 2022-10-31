function plotStatsLocs_upRegs(propIn,varIn,minArea,varName,lonLims,latLims,figname,classTypes,colmapCC)

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');
colormap('jet');

% Loop through projects
projects=fields(propIn);

for ii=1:length(projects)

    for jj=1:length(classTypes)
        thisTable=propIn.(projects{ii}).(classTypes{jj});
        thisTable(thisTable.area<minArea,:)=[];
        thisVar=thisTable.(varIn);
        s=subplot(4,3,jj);

        hold on

        if ~all(isnan(thisVar))
           scatter(thisTable.lon,thisTable.lat,25,thisVar,'filled');

            if length(thisVar)>15
                perc=prctile(thisVar,[5,95]);
                if length(unique(perc))>1
                    caxis(perc);
                end
            end
            colorbar
        end

        xlim(lonLims(ii,:));
        ylim(latLims(ii,:));

        plot(coastlon,coastlat,'-k')
        title([classTypes{jj}],'Color',colmapCC(jj,:));

        grid on
        box on
    end

    sgtitle([(projects{ii}),': ',varName],'FontSize',14,'FontWeight','bold');
    
    set(gcf,'PaperPositionMode','auto')
    print([figname,(projects{ii}),'.png'],'-dpng','-r0');
end
end