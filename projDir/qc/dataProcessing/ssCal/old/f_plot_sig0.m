function f_plot_sig0(PLT,outfile,xlims,titleIn,sfcWdspd,sfcWddir,avg_hdg, avg_alt,dB_bias_lb,attLiebe)
%plot sig0 of ocean scan data
%divide into rotation angles bigger or smaller than 180
sig0smaller=PLT.sig0(find(PLT.rota<=180));
elevSmaller=PLT.elev(find(PLT.rota<=180));

sig0bigger=PLT.sig0(find(PLT.rota>180));
elevBigger=PLT.elev(find(PLT.rota>180));

f2 = figure;
set(gcf,'Position',[200 500 800 600]);
%plot(PLT.elev,sig0);
hold on
plot(elevSmaller,sig0smaller,'r');
plot(elevBigger,sig0bigger,'b');
set(gca,'XLim',xlims);
set(gca,'YLim',[-20,20]);
xlabel('Incidence Angle, degrees');
ylabel('Normalized Radar Cross Section, dB');
title(titleIn);
set(0,'DefaultTextInterpreter','none');
tstr = sprintf(['Sfc wind sp = %s ms^-1\nSfc wind dir = %s deg\nAvg Hdg  = %5.0f deg\nAvg Alt = ' ...
    '%7.0f m\nLiebe bias = %4.2f dB\n2way attn = %5.2f dB'], num2str(round(sfcWdspd)), num2str(round(sfcWddir)),...
    avg_hdg, avg_alt,dB_bias_lb, 2*attLiebe);
text( xlims(2)-6,14,tstr);
legend('Rot angle <= 180 deg','Rot angle > 180 deg','location','southwest');
set(gcf,'PaperPositionMode','auto')
print(outfile,'-dpng','-r0')
end

