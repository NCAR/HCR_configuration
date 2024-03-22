% Width test

clear all;
close all;

sd1=5;
sd2=8;
mu=0;

x=-25:0.05:25;

y1=1/(sd1*sqrt(2*pi)).*exp(-0.5.*(x-mu).^2./sd1^2);
y2=1/(sd2*sqrt(2*pi)).*exp(-0.5.*(x-mu).^2./sd2^2);
y3=y1+y2;
y4=conv(y1,y3,'same');

plot(x,y1);
hold on
plot(x,y3);
plot(x,y4);

legend('y1','y3','y4')

% VEL
vel1=sum(y1.*x)/sum(y1);

vel3=sum(y3.*x)/sum(y3);
vel4=sum(y2.*x)/sum(y4);

% WIDTH
width1=(sum(y1.*(x-vel1).^2)./sum(y1)).^0.5;

width3=(sum(y3.*(x-vel3).^2)./sum(y3)).^0.5;
width4=(sum(y2.*(x-vel2).^2)./sum(y2)).^0.5;