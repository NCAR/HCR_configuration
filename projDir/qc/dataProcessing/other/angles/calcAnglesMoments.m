% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath('/h/eol/romatsch/codes/matlab/hcrQC/functions/');

%infile='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/20180116/cfrad.20180116_034400.165_to_20180116_034500.078_HCR_SUR.nc';
%infile='/scr/rain1/rsfdata/projects/socrates/hcr/qc/data/socrates/temp/20180116/cfrad.20180116_034400.165_to_20180116_034500.078_HCR_SUR.nc';
infile='/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/moments/10hz/20190807/cfrad.20190807_155000.100_to_20190807_155500.000_HCR_otrec.nc';

close all;

roll=ncread(infile,'roll');
pitch=ncread(infile,'pitch');
heading=ncread(infile,'heading');
rotation=ncread(infile,'rotation');
tilt=ncread(infile,'tilt');
drift=ncread(infile,'drift');

az_orig=ncread(infile,'azimuth');
el_orig=ncread(infile,'elevation');

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
%ylim([-91 -89])
legend('Orig')
title('Elevation angle [deg]');
subplot(4,1,2)
plot(azEl_cfrad(:,2))
%ylim([-91 -89])
legend('cfrad')
subplot(4,1,3)
plot(azEl_lee(:,2))
%ylim([-91 -89])
legend('Lee')
subplot(4,1,4)
hold on
plot(azEl_field(:,2))
yyaxis right
plot(roll)
%ylim([-91 -89])
legend('field')

disp(['Diff cfrad - Lee azimuth: ', num2str(max(azEl_cfrad(:,1)-azEl_lee(:,1)))]);
disp(['Diff cfrad - Lee elevation: ', num2str(max(azEl_cfrad(:,2)-azEl_lee(:,2)))]);

disp(['Diff cfrad - field azimuth: ', num2str(max(azEl_cfrad(:,1)-azEl_field(:,1)))]);
disp(['Diff cfrad - field elevation: ', num2str(max(azEl_cfrad(:,2)-azEl_field(:,2)))]);

disp(['Diff cfrad - orig azimuth: ', num2str(max(azEl_cfrad(:,1)-az_orig))]);
disp(['Diff cfrad - orig elevation: ', num2str(max(azEl_cfrad(:,2)-el_orig))]);