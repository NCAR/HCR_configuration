function f_plot_sig0measured(PLT,outfile,xlims,ylims,titleIn,wdSpeed,SST,windDir,hdgMinus180,atten,good)
%plot sig0 of ocean scan data
%divide into rotation angles bigger or smaller than 180
sig0smaller=PLT.sig0measured(find(PLT.rota<=180));
elevSmaller=PLT.elev(find(PLT.rota<=180));

sig0bigger=PLT.sig0measured(find(PLT.rota>180));
elevBigger=PLT.elev(find(PLT.rota>180));

f2 = figure;
set(gcf,'Position',[200 500 800 600]);
%plot(PLT.elev,sig0);
hold on
if good==1
    l1=plot(elevSmaller,sig0smaller,'r');
    l2=plot(elevBigger,sig0bigger,'b');
elseif good==0
    l1=plot(elevSmaller,sig0smaller,'m');
    l2=plot(elevBigger,sig0bigger,'color',[0.5 0.5 0.5]);
end
if isfield(PLT,'sig0model') & max(~isnan(PLT.sig0model))>0
    %plot(PLT.elev,PLT.sig0model(:,1),'c');
    l3=plot(PLT.elev,PLT.sig0model(:,2),'c','linewidth',2);
    %plot(PLT.elev,PLT.sig0model(:,3),'c');
    %plot(PLT.elev,PLT.sig0model(:,4),'g');
    l4=plot(PLT.elev,PLT.sig0model(:,5),'color',[0 0.5 0],'linewidth',2);
    %plot(PLT.elev,PLT.sig0model(:,6),'g');
    %plot(PLT.elev,PLT.sig0model(:,7),'k');
    l5=plot(PLT.elev,PLT.sig0model(:,8),'g','linewidth',2);
    %plot(PLT.elev,PLT.sig0model(:,9),'k');
    %l6=plot(PLT.polyElev,PLT.polySig0Model,'g','linewidth',3);
    plot(PLT.polyElev,PLT.polySig0Meas,'color',[1 1 1],'linewidth',4);
    l7=plot(PLT.polyElev,PLT.polySig0Meas,'color',[0 0 0],'linewidth',2);
    if ~isnan(SST)
        text(xlims(2)-10,ylims(2)-1.5,['Sea surf temp: ',num2str(SST),' C'],'fontsize',12);
    end
end
if ~isnan(wdSpeed)
    text(xlims(2)-10,ylims(2)-3,['Wind speed: ',num2str(wdSpeed),' m s^{-1}'],'fontsize',12);
end
if ~isnan(atten)
    text(xlims(2)-10,ylims(2)-4.5,['Two way atten: ',num2str(atten*2),' dB'],'fontsize',12);
end
if ~isnan(windDir)
    text(xlims(2)-10,ylims(1)+3,['Wind direction: ',num2str(windDir),' deg'],'fontsize',12);
end
if ~isnan(hdgMinus180)
    text(xlims(2)-10,ylims(1)+1.5,['Heading-180: ',num2str(hdgMinus180),' deg'],'fontsize',12);
end
set(gca,'XLim',xlims);
set(gca,'YLim',ylims);
xlabel('Elevation angle [deg]');
ylabel('Radar cross section [dB]');
title(titleIn,'interpreter','none');
if isfield(PLT,'sig0model') & max(~isnan(PLT.sig0model))>0
    %legend([l1,l2,l3,l4,l5],{'Rot angle <= 180 deg','Rot angle > 180 deg','FreilichVanhoff','Wu','CoxMunk'},'location','southwest');
    %legend([l1,l2,l6,l7],{'Rot angle <= 180 deg','Rot angle > 180 deg','CoxMunkFit','dataFit'},'location','southwest');
    legend([l1,l2,l3,l4,l5,l7],{'Rot angle <= 180 deg','Rot angle > 180 deg','FreilichVanhoff','Wu','CoxMunk','dataFit'},'location','southwest');
else
    legend([l1,l2],{'Rot angle <= 180 deg','Rot angle > 180 deg'},'location','southwest');
end
grid on
set(gcf,'PaperPositionMode','auto')
print(outfile,'-dpng','-r0')
end

