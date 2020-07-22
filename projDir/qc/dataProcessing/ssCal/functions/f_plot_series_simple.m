function f_plot_series_simple(PLT,outfile)
%plot time series of ocean scan data
h = [];
f = figure;
set(gcf,'Position',[1 100 1600 900 ]);

h(1) = subplot(4,1,1);
hold on;
plot(PLT.time,PLT.rota);
xlim([PLT.time(1) PLT.time(end)]);
title('Rotation angle and surface velocity');
ylabel('[deg]');

yyaxis right;
plot(PLT.time,PLT.vel);
ylabel('Vel [m/s]');
hold off;

h(2) = subplot(4,1,2);
plot(PLT.time,PLT.refl);
xlim([PLT.time(1) PLT.time(end)]);
title('Maximum reflectivity');
ylabel('[dBZ]');

h(3) = subplot(4,1,3);
plot(PLT.time,PLT.pitch);
hold on;
plot(PLT.time,PLT.roll);
title('Pitch and roll angles');
ylabel('[deg]');
legnd = [{'Pitch'}; {'Roll'}];
xlim([PLT.time(1) PLT.time(end)]);
legend( legnd, 'Location', 'NorthEast', 'Interpreter', 'none' );

h(4) = subplot(4,1,4);
plot(PLT.time,PLT.alt);
hold on;
plot(PLT.time,PLT.range);
hold off;
linkaxes(h,'x');
title('Aircraft altitude and range to max dBZ');
ylabel('[m]');
legnd = [ {'Aircraft'}; {'Max dBZ'} ];
legend( legnd, 'Location', 'NorthEast', 'Interpreter', 'none' );
xlim([PLT.time(1) PLT.time(end)]);

set(gcf,'PaperPositionMode','auto');
print(outfile,'-dpng','-r0');
end

