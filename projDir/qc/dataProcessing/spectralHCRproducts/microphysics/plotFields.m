function plotFields(dataIn,fieldsIn,figdir,ylims,showPlot)

clims.DBZ=[-40,30];
clims.LDR=[-35,-5];
clims.VEL=[-15,15];
clims.WIDTH=[0,3];
clims.SKEW=[-3,3];
clims.KURT=[-6,6];
clims.EE_WIDTH=[0,13];
clims.L_SLOPE=[0,20];
clims.R_SLOPE=[-20,0];
clims.MELTING_LAYER=[6,24];

col1=cat(1,[0,0,0],jet);
col2=cat(1,[0,0,0],velCols);

col.DBZ=col1;
col.LDR=col1;
col.VEL=col2;
col.WIDTH=col1;
col.SKEW=col2;
col.KURT=col2;
col.EE_WIDTH=col1;
col.L_SLOPE=col1;
col.R_SLOPE=flipud(col1);
col.MELTING_LAYER=flipud(col1);

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
print(f1,[figdir,'inFields_',datestr(dataIn.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(dataIn.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end