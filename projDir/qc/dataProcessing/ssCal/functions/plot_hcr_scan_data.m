
% just the ploting utility that goes with determine_bias and
% plot_sfc_info.  Note that special consideration is given to 
% allow rational re-labeling of the time axis during window
% resizing (this resizing requires a global variable and the 
% myXtickResize and label_time_axis subroutines).

h = [];
   f = figure;
set(gcf,'Position',[1 100 1600 900 ]);

   h(1) = subplot(4,1,1);
   hold on;
   yyaxis left;
   plot(PLT.time,PLT.rota);
   tlims = get(gca,'XLim');
   %   label_time_axis();
   title('Rotation angle and surface velocity');
   ylabel('[deg]');

   yyaxis right;
   plot(PLT.time,PLT.vel);
   ylabel('Vel [m/s]');
   hold off;

   h(2) = subplot(4,1,2);
   plot(PLT.time,PLT.refl);
   set(gca,'XLim',tlims); % reflectivity can have a lot of NaNs
   %   label_time_axis();
   title('Maximum reflectivity');
   ylabel('[dBZ]');

   h(3) = subplot(4,1,3);
   plot(PLT.time,PLT.pitch);
   hold on;
   plot(PLT.time,PLT.roll);
   plot(PLT.time,PLT.angdif);
   %   label_time_axis();
   title('Pitch and roll angles, elevation angle difference');
   ylabel('[deg]');
   legnd = [{'Pitch'}; {'Roll'}; {'Rot+roll+pitch - elev'}];
   legend( legnd, 'Location', 'NorthEast', 'Interpreter', 'none' );

   h(4) = subplot(4,1,4);
   plot(PLT.time,PLT.alt);
   hold on;
   plot(PLT.time,PLT.range);
   hold off;
   linkaxes(h,'x');
   %label_time_axis();
   title('Aircraft altitude and range to max dBZ');
   ylabel('[m]');
   legnd = [ {'Aircraft'}; {'Max dBZ'} ];
   legend( legnd, 'Location', 'NorthEast', 'Interpreter', 'none' );


