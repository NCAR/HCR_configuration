% Membership coefficients

dbz.rain=[1,3]; % Rain threshold of 5 by shupe 2007
dbz.drizzle=[-15,-13,5,7]; % Drizzle, cloud threshold of -15 by  Measurement of Stratus Cloud and Drizzle Parameters in ASTEX with a Kα-Band Doppler Radar and a Microwave Radiometer A. S. Frisch1, C. W. Fairall1, and J. B. Snider1
dbz.cloud=[-11,-9];
dbz.mixed=[-3,-1];
dbz.lfrozen=[2,4]; % Snow threshold of 5 by shupe 2007
dbz.sfrozen=[6,8];

ldr.rain=[-26,-24];
ldr.drizzle=[-24,-22];
ldr.cloud=[-24,-22];
ldr.mixed=[-20,-17,-8,-6];
ldr.lfrozen=[-29,-27,-15,-12];
ldr.sfrozen=[-27,-25,-15,-12];

vel.rain=[3,3.5];
vel.drizzle=[0,0.5,3.5,4];
vel.cloud=[0.5,1];
vel.mixed=[0.5,1];
vel.lfrozen=[0.5,1];
vel.sfrozen=[1,2];

temp.rain=[-80,-10];
temp.drizzle=[-90,-10];
temp.cloud=[-100,-10];
temp.mixed=[-2,0,3,6];
temp.lfrozen=[0,6];
temp.sfrozen=[0,5];