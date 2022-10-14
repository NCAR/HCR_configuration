function plotLocsProjects(lons,lats,lonLims,latLims,figname,classTypes,colmapCC)

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,1200],'renderer','painters','visible','on');

% CSET
s1=subplot(3,1,1);
hold on
for jj=1:length(classTypes)
    thisLons=lons.cset.(classTypes{jj});
    thisLats=lats.cset.(classTypes{jj});

    scatter(thisLons,thisLats,'+','MarkerEdgeColor',colmapCC(jj,:));
end

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')

grid on
box on

title('CSET')

s1pos=s1.Position;

ldg=legend(classTypes,'Location','southoutside');
ldg.NumColumns=4;

s1.Position=s1pos;

% SOCRATES

s2=subplot(3,1,2);
hold on
for jj=1:length(classTypes)
    thisLons=lons.socrates.(classTypes{jj});
    thisLats=lats.socrates.(classTypes{jj});

    scatter(thisLons,thisLats,'+','MarkerEdgeColor',colmapCC(jj,:));
end

xlim(lonLims(2,:));
ylim(latLims(2,:));

plot(coastlon,coastlat,'-k')

grid on
box on

title('SOCRATES')

% OTREC

s3=subplot(3,1,3);
hold on
for jj=1:length(classTypes)
    thisLons=lons.otrec.(classTypes{jj});
    thisLats=lats.otrec.(classTypes{jj});

    scatter(thisLons,thisLats,'+','MarkerEdgeColor',colmapCC(jj,:));
end

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')

grid on
box on

title('OTREC')

set(gcf,'PaperPositionMode','auto')
print(figname,'-dpng','-r0');
end