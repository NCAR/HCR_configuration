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
    figs.(['fig',num2str(ii)]).PaperPositionMode = 'manual';
    figs.(['fig',num2str(ii)]).PaperUnits = 'inches';
    figs.(['fig',num2str(ii)]).Units = 'inches';
    figs.(['fig',num2str(ii)]).PaperPosition = [0, 0, wi, hi];
    figs.(['fig',num2str(ii)]).PaperSize = [wi, hi];
    figs.(['fig',num2str(ii)]).Resize = 'off';
    figs.(['fig',num2str(ii)]).InvertHardcopy = 'off';
    
    set(figs.(['fig',num2str(ii)]),'color','w');
    
    for jj=1:numSub(ii)
        ax1=subplot(numSub(ii),1,jj,'units','inch');
        posInd=numSub(ii)-jj;
        fig1=surf(data.time,asl,data.(fields{kk}));
        ylim(yLimits);
        ax1.Position = [0.5+1/(yLimits(2)-yLimits(1))/2,hi/numSub(ii)*posInd+0.5,...
            wi-0.6-1/(yLimits(2)-yLimits(1))/1.5,hi/numSub(ii)-0.8]; % We have to adjust for yLimits because that changes the plot width
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

