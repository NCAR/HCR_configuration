function f_plot_series(PLT,outfile)
%plot time series of ocean scan data
plot_hcr_scan_data;
set(gcf,'PaperPositionMode','auto');
print(outfile,'-dpng','-r0');
end

