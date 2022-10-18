function plotLocsTypes(lons,lats,lonLims,latLims,figname,classTypes,colmapCC)

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');

% Loop through projects
projects=fields(lons);

for ii=1:length(projects)

    for jj=1:length(classTypes)
        s=subplot(4,3,jj);

        hold on

        thisLons=lons.(projects{ii}).(classTypes{jj});
        thisLats=lats.(projects{ii}).(classTypes{jj});

        scatter(thisLons,thisLats,'+','MarkerEdgeColor',colmapCC(jj,:));

        xlim(lonLims(ii,:));
        ylim(latLims(ii,:));

        plot(coastlon,coastlat,'-k')
        title([classTypes{jj}],'Color',colmapCC(jj,:));

        grid on
        box on
    end

    sgtitle((projects{ii}),'FontSize',14,'FontWeight','bold');
    
    set(gcf,'PaperPositionMode','auto')
    print([figname,(projects{ii}),'.png'],'-dpng','-r0');
end
end