function [fig,s]=do_plotHCR(data,ylimUpper);
% Plot HCR variables

infields=fields(data);
plotVars=length(infields)-7;

% Get indices
numInds=floor(size(data.(infields{1}),2)./5000);
plotInds=1:max([numInds,1]):size(data.(infields{1}),2);

fig=figure('Position',[200 500 1200 min([plotVars*350,1700])],'DefaultAxesFontSize',12);
colormap('jet');

for ii=1:plotVars
    s.(infields{ii})=subplot(plotVars,1,ii);
    surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.(infields{ii})(:,plotInds),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    box on
    title(infields{ii},'interpreter','none');
end
end