% Membership coefficients

dbz.rain=[1,3]; % Rain threshold of 5 by shupe 2007
dbz.drizzle=[-15,-13,5,7]; % Drizzle, cloud threshold of -15 by  Measurement of Stratus Cloud and Drizzle Parameters in ASTEX with a Kα-Band Doppler Radar and a Microwave Radiometer A. S. Frisch1, C. W. Fairall1, and J. B. Snider1
dbz.cloud=[-11,-9];
dbz.mixed=[];
dbz.lfrozen=[2,4]; % Snow threshold of 5 by shupe 2007
dbz.sfrozen=[6,8];

hcrldr.rain=[-26,-24];
hcrldr.drizzle=[-24,-22];
hcrldr.cloud=[-24,-22];
hcrldr.mixed=[-16,-14,-8,-6];
hcrldr.lfrozen=[-29,-27,-14,-13];
hcrldr.sfrozen=[-28,-26,-14,-13];

vel.rain=[3,3.5];
vel.drizzle=[0,0.5,3.5,4];
vel.cloud=[0.5,1];
vel.mixed=[-1.5,-1];
vel.lfrozen=[0.5,1];
vel.sfrozen=[1,2];

temp.rain=[-80,-10];
temp.drizzle=[-90,-10];
temp.cloud=[-100,-10];
temp.mixed=[-8,-6,3,6];
temp.lfrozen=[0,6];
temp.sfrozen=[0,5];

back.rain=[1e-6,1e-5];
back.drizzle=[1e-6,1e-5];
back.cloud=[1e-7,1e-6];
back.mixed=[];
back.lfrozen=[1e-5,1e-4];
back.sfrozen=[1e-6,1e-5];

hsrlldr.rain=[0.1,0.2];
hsrlldr.drizzle=[0.1,0.2];
hsrlldr.cloud=[0.15,0.25];
hsrlldr.mixed=[0.1,0.2];
hsrlldr.lfrozen=[0,0.1];
hsrlldr.sfrozen=[0.05,0.15];