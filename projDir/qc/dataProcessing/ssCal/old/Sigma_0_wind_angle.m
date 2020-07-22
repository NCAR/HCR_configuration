close all; clear; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lihua Li Seattle AMS Radar conference  pp 204 -207, 6-12, August 2003
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref_index= complex(7.56, -13.6);  % Refractive index for 20 C sea surafce temperature

eff_ref_coeff= 0.89*(ref_index-1)/(ref_index+1);  % 0.89 Ocean surface roughness correction factor Freilich, 2003

fac1=abs(eff_ref_coeff)^2;


scan_angle=[-30:0.2:30]*pi/180;
surface_wind=[0:0.05:15]; % m/sec
[angle, wind]=meshgrid(scan_angle,surface_wind);
sur_slope= (0.003+1.92e-3*wind)+ (3.16e-3*wind); % Cox and Munk 1954

fac2= sur_slope.*(cos(angle)).^4;

fac3=exp(-(tan(angle)).^2./sur_slope);


sigma_0=10.*log10(((fac1*fac3./fac2)));

figure; contourf(angle*180/pi, wind,sigma_0); colorbar; 
xlabel('Angles')
ylabel('Surface wind  m s^{-1}')



