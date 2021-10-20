% Membership coefficients

dbz.rain=[2,4]; % Rain threshold of 5 by shupe 2007
dbz.drizzle=[-14,-12,4,6]; % Drizzle, cloud threshold of -15 by  Measurement of Stratus Cloud and Drizzle Parameters in ASTEX with a Kα-Band Doppler Radar and a Microwave Radiometer A. S. Frisch1, C. W. Fairall1, and J. B. Snider1
dbz.cloud=[-12,-10,];
dbz.mixed=[-3,-1];
dbz.lfrozen=[3,5]; % Snow threshold of 5 by shupe 2007
dbz.sfrozen=[5,7];

ldr.rain=[-26,-24];
ldr.drizzle=[-25,-23];
ldr.cloud=[-25,-23];
ldr.mixed=[-20,-17,-8,-6];
ldr.lfrozen=[-28,-26,-15,-12];
ldr.sfrozen=[-27,-25,-15,-12];

vel.rain=[3,3.5];
vel.drizzle=[0,0.5,3.5,4];
vel.cloud=[0.5,1];
vel.mixed=[0.5,1];
vel.lfrozen=[0.5,1];
vel.sfrozen=[-1,0,1,2];

width.rain=[0.2,0.4];
width.drizzle=[0.2,0.4];
width.cloud=[0.2,0.4];
width.mixed=[0.5,0.6];
width.lfrozen=[0.4,0.6];
width.sfrozen=[0.4,0.6];

temp.rain=[-80,-10];
temp.drizzle=[];
temp.cloud=[];
temp.mixed=[-2,0,3,6];
temp.lfrozen=[0,6];
temp.sfrozen=[0,5];