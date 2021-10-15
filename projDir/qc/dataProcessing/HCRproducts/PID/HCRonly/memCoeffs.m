% Membership coefficients

dbz.rain=[3,5]; % Rain threshold of 5 by shupe 2007
dbz.drizzle=[-17,-15,5,8]; % Drizzle, cloud threshold of -15 by  Measurement of Stratus Cloud and Drizzle Parameters in ASTEX with a Kα-Band Doppler Radar and a Microwave Radiometer A. S. Frisch1, C. W. Fairall1, and J. B. Snider1
dbz.cloud=[-15,-13,];
dbz.mixed=[-3,-1,20,25];
dbz.lfrozen=[3,5,18,20]; % Snow threshold of 5 by shupe 2007
dbz.sfrozen=[5,7];

ldr.rain=[-27,-22];
ldr.drizzle=[-27,-25];
ldr.cloud=[-27,-25];
ldr.mixed=[-20,-17,-8,-6];
ldr.lfrozen=[-22,-20,-16,-14];
ldr.sfrozen=[-26,-24,-15,-12];

vel.rain=[2.5,3];
vel.drizzle=[0.5,1,3,3.5];
vel.cloud=[1,1.5];
vel.mixed=[0.5,1,3,4];
vel.lfrozen=[0.5,1];
vel.sfrozen=[-1,0,1,2];

width.rain=[0.1,0.2];
width.drizzle=[0.2,0.3];
width.cloud=[0.1,0.2];
width.mixed=[0.2,0.3];
width.lfrozen=[0.2,0.3];
width.sfrozen=[0.7,0.9];

temp.rain=[-2,2];
temp.drizzle=[-50,-39];
temp.cloud=[];
temp.mixed=[-2,0,3,6];
temp.lfrozen=[0,6];
temp.sfrozen=[-1,5];