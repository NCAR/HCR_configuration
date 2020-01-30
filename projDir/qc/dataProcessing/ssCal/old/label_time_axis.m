function label_time_axis();

% works only in conjunction with the custom routine myXtickResize

% Re-labels the X-axis of a plot using string equivalents for 
% the time; assumes data are plotted using time, in seconds.
% The time axis should be relabeled whenever the plot is created,
% zoomed, or resized, and this routine takes care of that.

% This routine will not work correctly if the "reset view" option is
% selected from a figure's Tools menu (not my fault -- Matlab's!)

% This routine may be used on multiple figure windows, without conflict.
% This is because any create/zoom/resize is generally a manual process,
% and the current figure (gcf, and gca) is set to that window during a 
% manual operation.

myXtickResize();

z = zoom(gcf);
set(z, 'ActionPostCallback',@myXtickResize);
z = zoom(gca);
set(z, 'ActionPostCallback',@myXtickResize);
set(gcf,'SizeChangedFcn',@myXtickResize);
set(gcf,'ResizeFcn',@myXtickResize);
set(gcf,'CreateFcn',@myXtickResize);
%set(gca,'CreateFcn',@myXtickResize);

pp = pan;
set(pp,'ActionPostCallback',@myXtickResize);

