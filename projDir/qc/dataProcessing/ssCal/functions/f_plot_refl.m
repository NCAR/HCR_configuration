function f_plot_refl(PLT,outfile,xlims,ylims,titleIn)
%plot sig0 of ocean scan data
%divide into rotation angles bigger or smaller than 180
reflsmaller=PLT.refl(find(PLT.rota<=180));
elevSmaller=PLT.elev(find(PLT.rota<=180));

reflbigger=PLT.refl(find(PLT.rota>180));
elevBigger=PLT.elev(find(PLT.rota>180));

f2 = figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(elevSmaller,reflsmaller,'r');
plot(elevBigger,reflbigger,'b');
set(gca,'XLim',xlims);
set(gca,'YLim',ylims);
xlabel('Elevation angle [deg]');
ylabel('Reflectivity [dBZ]');
title(titleIn,'interpreter','none');
legend('Rot angle <= 180 deg','Rot angle > 180 deg','location','southwest');
grid on
set(gcf,'PaperPositionMode','auto')
print(outfile,'-dpng','-r0')
end

