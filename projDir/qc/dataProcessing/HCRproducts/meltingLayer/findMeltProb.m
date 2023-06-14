function meltProb=findMeltProb(data,velDiff,dbzDiff)
% Fuzzy logic scheme to determine the probability for melting layer
% TEMP, LDR, DBZ, velDiff, dbzDiff
w=[0.1,0.21,0.07,0.39,0.23];

m(1,:,:)=trapmf(data.TEMP,[-1,0,1,7]);
m(2,:,:)=smf(data.LDR,[-20,-18]);
m(3,:,:)=smf(data.DBZ_MASKED,[-20,10]);
m(4,:,:)=smf(velDiff,[0.1,0.4]);
m(5,:,:)=trapmf(dbzDiff,[0.3,0.5,3,4]);

mw=m.*w';
meltProb=squeeze(sum(mw,1));
end