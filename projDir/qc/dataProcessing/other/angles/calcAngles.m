% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath('/h/eol/romatsch/gitPriv/process_HCR/hcrQC/functions/');

infile='/scr/rain1/rsfdata/projects/socrates/hcr/qc/data/socrates/SOCRATES.rf01.20180115_215100.to.20180116_052600.txt';

[header indata]=robustcsvreadheader(infile,' ');

header(1)=[];

header_new = header;
uniqStr = unique(header);
for ii = 1:length(uniqStr)
    idx = strmatch(uniqStr{ii}, header, 'exact');
    for jj = 2:length(idx)
        header_new{idx(jj)} = [header{idx(jj)}, '_', num2str(jj)];
    end
end

intable=array2table(indata,'VariableNames',header_new);

%% Calculate
close all;

roll=intable.rollSec;
pitch=intable.pitchSec;
heading=intable.headingSec;
rotation=intable.rotation;
tilt=intable.tilt;
drift=intable.drift;

az_orig=intable.azimuth;
el_orig=intable.elevation;

[azEl_cfrad,azEl_lee,azEl_field] = azEl(roll,pitch,heading,rotation,tilt,drift);

hFig=figure;
set(hFig, 'Position', [100 100 2800 1200])
subplot(4,1,1)
plot(az_orig);
legend('Orig')
title('Azimuth [deg]');
subplot(4,1,2)
plot(azEl_cfrad(:,1))
legend('cfrad')
subplot(4,1,3)
plot(azEl_lee(:,1))
legend('Lee')
subplot(4,1,4)
plot(azEl_field(:,1))
legend('field')

hFig2=figure;
set(hFig2, 'Position', [100 100 2800 1200])
subplot(4,1,1)
plot(el_orig);
ylim([-91 -89])
legend('Orig')
title('Elevation angle [deg]');
subplot(4,1,2)
plot(azEl_cfrad(:,2))
ylim([-91 -89])
legend('cfrad')
subplot(4,1,3)
plot(azEl_lee(:,2))
ylim([-91 -89])
legend('Lee')
subplot(4,1,4)
plot(azEl_field(:,2))
ylim([-91 -89])
legend('field')

disp(['Diff cfrad - Lee azimuth: ', num2str(max(azEl_cfrad(:,1)-azEl_lee(:,1)))]);
disp(['Diff cfrad - Lee elevation: ', num2str(max(azEl_cfrad(:,2)-azEl_lee(:,2)))]);

disp(['Diff cfrad - field azimuth: ', num2str(max(azEl_cfrad(:,1)-azEl_field(:,1)))]);
disp(['Diff cfrad - field elevation: ', num2str(max(azEl_cfrad(:,2)-azEl_field(:,2)))]);

disp(['Diff cfrad - orig azimuth: ', num2str(max(azEl_cfrad(:,1)-az_orig))]);
disp(['Diff cfrad - orig elevation: ', num2str(max(azEl_cfrad(:,2)-el_orig))]);