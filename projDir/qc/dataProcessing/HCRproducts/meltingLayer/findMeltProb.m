function meltProb=findMeltProb(data,velDiff)
% Fuzzy logic scheme to determine the probability for melting layer
% DBZ, TEMP, LDR, velDiff
w=[0.15,0.25,0.25,0.35];

m(1,:,:)=smf(data.DBZ_MASKED,[-30,10]);
m(2,:,:)=trapmf(data.TEMP,[-1,0,1,7]);
m(3,:,:)=smf(data.LDR,[-20,-18]);
m(4,:,:)=smf(velDiff,[0.1,0.4]);

mw=m.*w';
meltProb=squeeze(sum(mw,1));
end