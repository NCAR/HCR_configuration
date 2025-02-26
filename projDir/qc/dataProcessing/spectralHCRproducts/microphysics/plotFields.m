function plotFields(dataIn,fieldsIn,type,campaign,figdir,ylims,showPlot)

clims.DBZ_MASKED=[-40,30];
clims.LDR=[-35,-5];
clims.VEL_MASKED=[-15,15];
clims.WIDTH_SPEC=[0,3];
clims.SKEWNESS=[-3,3];
clims.KURTOSIS=[-6,6];
clims.EDGE_EDGE_WIDTH=[0,13];
clims.LEFT_SLOPE=[0,20];
clims.RIGHT_SLOPE=[-20,0];
clims.MELTING_LAYER=[6,24];
clims.TEMP=[-50,30];

col1=cat(1,[0,0,0],jet);
col2=cat(1,[0,0,0],velCols);
col3=flipud(cat(1,jet,[0,0,0]));

col.DBZ_MASKED=col1;
col.LDR=col1;
col.VEL_MASKED=col2;
col.WIDTH_SPEC=col1;
col.SKEWNESS=col2;
col.KURTOSIS=col2;
col.EDGE_EDGE_WIDTH=col1;
col.LEFT_SLOPE=col1;
col.RIGHT_SLOPE=col3;
col.MELTING_LAYER=col3;
col.TEMP=col1;

rows=min([4,ceil(sqrt(length(fieldsIn)))]);
cols=ceil(length(fieldsIn)/rows);

%% Figure
f1 = figure('Position',[200 500 cols*700 rows*300],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(rows,cols,'TileSpacing','tight','Padding','tight');

for ii=1:length(fieldsIn)

    s1=nexttile(ii);

    dataIn.(fieldsIn{ii})(isnan(dataIn.(fieldsIn{ii})))=-999;

    hold on
    surf(dataIn.time,dataIn.asl./1000,dataIn.(fieldsIn{ii}),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim(clims.(fieldsIn{ii}));
    s1.Colormap=col.(fieldsIn{ii});
    colorbar
    grid on
    box on
    title(fieldsIn{ii},'interpreter','none')
    ylim(ylims);
    xlim([dataIn.time(1),dataIn.time(end)]);

end

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'inFields_',type,'_',campaign,'_', ...
    datestr(dataIn.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(dataIn.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end