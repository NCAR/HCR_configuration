function plotLocsTypes(lons,lats,lonLims,latLims,figname,classTypes,colmapCC)

load coastlines

%% Scatter
fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');

% Loop through projects
projects=fields(lons);

for bb=1:length(projects)

    for aa=1:length(classTypes)
        s=subplot(4,3,aa);

        hold on

        thisLons=lons.(projects{bb}).(classTypes{aa});
        thisLats=lats.(projects{bb}).(classTypes{aa});

        scatter(thisLons,thisLats,'+','MarkerEdgeColor',colmapCC(aa,:));

        xlim(lonLims(bb,:));
        ylim(latLims(bb,:));

        plot(coastlon,coastlat,'-k')
        title([classTypes{aa}],'Color',colmapCC(aa,:));

        grid on
        box on
    end

    sgtitle((projects{bb}),'FontSize',14,'FontWeight','bold');
    
    set(gcf,'PaperPositionMode','auto')
    print([figname,(projects{bb}),'.png'],'-dpng','-r0');
end

%% Per flight hour
load('/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/flightHourGrids.mat');

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1800,1200],'renderer','painters','visible','on');
colormap('jet');

% Loop through projects
projects=fields(lons);
gridStep=1;

for bb=1:length(projects)

    lonLength=(lonLims(bb,2)-lonLims(bb,1))/gridStep;
    latLength=(latLims(bb,2)-latLims(bb,1))/gridStep;
    
    lonSteps=lonLims(bb,1):gridStep:lonLims(bb,2);
    latSteps=latLims(bb,1):gridStep:latLims(bb,2);

    for aa=1:length(classTypes)
        hourGrid=zeros(latLength,lonLength);
        s=subplot(4,3,aa);

        hold on

        thisLons=lons.(projects{bb}).(classTypes{aa});
        thisLats=lats.(projects{bb}).(classTypes{aa});

        % Loop through output grid
        for ii=1:size(hourGrid,1)
            for jj=1:size(hourGrid,2)
                pixInds=find(thisLats>latSteps(ii) & thisLats<=latSteps(ii+1) & ...
                    thisLons>lonSteps(jj) & thisLons<=lonSteps(jj+1));
                hzPix=length(pixInds);
                hourGrid(ii,jj)=hourGrid(ii,jj)+hzPix;
            end
        end

        perHour=hourGrid./flightHourGrids.(projects{bb});
        perHour(perHour==0)=nan;

        h=imagesc(lonSteps(1:end-1)+gridStep/2,latSteps(1:end-1)+gridStep/2,perHour);
        set(h,'AlphaData',~isnan(perHour));

        caxis([0,10]);
        colorbar;

        xlim(lonLims(bb,:));
        ylim(latLims(bb,:));

        plot(coastlon,coastlat,'-k')
        title([classTypes{aa}]);

        grid on
        box on
    end

    sgtitle((projects{bb}),'FontSize',14,'FontWeight','bold');
    
    set(gcf,'PaperPositionMode','auto')
    print([figname,(projects{bb}),'_perFlightHour.png'],'-dpng','-r0');
end
end