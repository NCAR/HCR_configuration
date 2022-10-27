function plotStatsLocs(varIn,varName,lons,lats,lonLims,latLims,figname,classTypes,colmapCC)

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');
colormap('jet');

% Loop through projects
projects=fields(lons);

for ii=1:length(projects)
   
    for jj=1:length(classTypes)
        s=subplot(4,3,jj);

        hold on

        thisLons=lons.(projects{ii}).(classTypes{jj});
        thisLats=lats.(projects{ii}).(classTypes{jj});

        scatter(thisLons,thisLats,25,varIn.(projects{ii}).(classTypes{jj}),'filled');

        if length(varIn.(projects{ii}).(classTypes{jj}))>15
            perc=prctile(varIn.(projects{ii}).(classTypes{jj}),[5,95]);
            if length(unique(perc))>1
                caxis(perc);
            end
        end
        colorbar

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