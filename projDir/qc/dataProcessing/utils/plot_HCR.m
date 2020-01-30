function figs = plot_HCR(data,fields,yLimits)

maxY=10;
plotsPerFig=3;

%Plot HCR variable
aslM=HCRrange2asl(data.range,data.elevation,data.altitude);
asl=aslM./1000;

numSub=repmat(plotsPerFig,1,floor(size(fields,1)/plotsPerFig));

remPlots=rem(size(fields,1),plotsPerFig);

if remPlots>0
    numSub=[numSub,remPlots];
end

kk=1;

for ii=1:length(numSub)
    wi=10;
    hi=numSub(ii)*4;
    
    figs.(['fig',num2str(ii)])=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[1,100,wi,hi]);
    
    for jj=1:numSub(ii)
        ax1=subplot(numSub(ii),1,jj);
        fig1=surf(data.time,asl,data.(fields{kk}));
        ylim(yLimits);
        fig1.EdgeColor='none';
        xlim([data.time(1),data.time(end)]);
        view(2);
        
        colormap(jet)
        hcb=colorbar;
        ylabel('Altitude [km]');
        title([fields{kk},': ',datestr(data.time(1),'YYYY-mm-DD hh:MM:ss'),' to ',...
            datestr(data.time(end),'YYYY-mm-DD hh:MM:ss')],'interpreter','none');
        
        kk=kk+1;
    end
    drawnow;
end
end

