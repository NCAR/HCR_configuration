function pix=hcrTimeToPix(timeMins,pixSec)
% Calculate the hcr pixel number that correspond to a timespan in minutes
% timeMins is the desired timespan
% pixSec is the time resolution of the data in seconds

pix=1/pixSec*60*timeMins;
end