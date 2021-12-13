function [fig,s]=do_plotHCR(data,ylimUpper);
% Plot HCR variables

infields=fields(data);
plotVars=length(infields)-7;

fig=figure('Position',[200 500 1200 plotVars*350],'DefaultAxesFontSize',12);
colormap('jet');

for ii=1:plotVars
    s.(infields{ii})=subplot(plotVars,1,ii);
    surf(data.time,data.asl./1000,data.(infields{ii}),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title(infields{ii});
end
end