function meltProb=findMeltProb(data,velDiff,dbzDiff)
% Fuzzy logic scheme to determine the probability for melting layer
% TEMP, LDR, DBZ, velDiff, dbzDiff
w=[0.1,0.21,0.07,0.39,0.23];

data.TEMP(isnan(data.TEMP))=-999;
data.LDR(isnan(data.LDR))=-999;
data.DBZ(isnan(data.DBZ))=-999;
velDiff(isnan(velDiff))=-999;
dbzDiff(isnan(dbzDiff))=-999;

m(1,:,:)=trapmf(data.TEMP,[-1,0,1,5]);
m(2,:,:)=trapmf(data.LDR,[-20,-18,-6,-5]);
m(3,:,:)=smf(data.DBZ,[-20,10]);
m(4,:,:)=smf(velDiff,[0.1,0.4]);
m(5,:,:)=trapmf(dbzDiff,[0.3,0.5,3,4]);

mw=m.*w';
meltProb=squeeze(sum(mw,1));

meltProb(data.TEMP==-999)=nan;
end