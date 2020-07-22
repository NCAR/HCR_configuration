
% just the ploting utility that goes with plot_sfc_info

h = [];
   f = figure;
set(gcf,'Position',[1 100 1600 900 ]);

   h(1) = subplot(4,1,1);
   plot(PLT.time,PLT.rota);
   tlims = get(gca,'XLim');
   %   label_time_axis();
   title('Rotation Angle');
   ylabel('degrees');

   h(2) = subplot(4,1,2);
   plot(PLT.time,PLT.refl);
   set(gca,'XLim',tlims); % reflectivity can have a lot of NaNs
   %   label_time_axis();
   title('Maximum Reflectivity');
   ylabel('dBZ');

   h(3) = subplot(4,1,3);
   plot(PLT.time,PLT.pitch);
   %   label_time_axis();
   title('Pitch Angle');
   ylabel('degrees');

   h(4) = subplot(4,1,4);
   plot(PLT.time,PLT.alt);
   hold on;
   plot(PLT.time,PLT.range);
   hold off;
   linkaxes(h,'x');
   label_time_axis();
   title('Altitude/Range');
   ylabel('meters');
   legnd = [ {'Aircraft'}; {'Max dBZ'} ];
   legend( legnd, 'Location', 'NorthEast', 'Interpreter', 'none' );

   subplot(h(1));
   tHandle = text( 'String', ['Data of ' dstr], 'Units', 'normalized', 'Position', [0.02,1.15]);

   % this statement ensures that Tools menu Reset View does work
   % (this
   % was an sob to figure out!):
gh = guihandles;
set (gh.figMenuResetView,'Callback',['toolsmenufcn ResetView; ' ...
                    'myXtickResize()']);

myXtickResize
