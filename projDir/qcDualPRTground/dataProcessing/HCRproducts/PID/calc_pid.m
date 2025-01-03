function[classOut]=calc_pid(data,plotIn)
%   Membership functions for particle detection

memCoeffs
% DBZ, LDR, VEL, TEMP
w=[34 22 22 22];

% pid_hcr (number before post processing)
%  1 Rain
%  2 Drizzle
%  3 Cloud liquid
%  4 Melting (is added in post processing)
%  5 Large frozen
%  6 Small frozen

dbzOrig=data.DBZ;
data.DBZ(isnan(data.DBZ))=-999;
data.VEL(isnan(data.VEL))=-999;
data.LDR(isnan(data.LDR))=-999;
data.TEMP(isnan(data.TEMP))=-999;

result=nan(5,size(data.LDR,1),size(data.LDR,2));

%  Membership functions for rain
m=nan(4,size(data.DBZ,1),size(data.DBZ,2));
m(1,:,:)=smf(data.DBZ,[dbz.rain(1),dbz.rain(2)]);  % Rain
m(2,:,:)=trapmf(data.LDR,[-100,-99,ldr.rain(1),ldr.rain(2)]);
m(3,:,:)=smf(data.VEL,[vel.rain(1),vel.rain(2)]); % Beard and Pruppacher 1969
m(4,:,:)=smf(data.TEMP,[temp.rain(1),temp.rain(2)]);

result(1,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(1,:,:),'Rain',plotIn);
end

%  Membership functions for drizzle
m=nan(4,size(data.DBZ,1),size(data.DBZ,2));
m(1,:,:)=trapmf(data.DBZ,[dbz.drizzle(1),dbz.drizzle(2),dbz.drizzle(3),dbz.drizzle(4)]);
m(2,:,:)=trapmf(data.LDR,[-100,-99,ldr.drizzle(1),ldr.drizzle(2)]);
m(3,:,:)=trapmf(data.VEL,[vel.drizzle(1),vel.drizzle(2),vel.drizzle(3),vel.drizzle(4)]);% Houze 2014, Beard and Pruppacher 1969
m(4,:,:)=smf(data.TEMP,[temp.drizzle(1),temp.drizzle(2)]);

result(2,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(2,:,:),'Drizzle',plotIn);
end

%  Membership functions for cloud liquid
m=nan(4,size(data.DBZ,1),size(data.DBZ,2));
m(1,:,:)=trapmf(data.DBZ,[-100,-99,dbz.cloud(1),dbz.cloud(2)]);
m(2,:,:)=trapmf(data.LDR,[-100,-99,ldr.cloud(1),ldr.cloud(2)]);
m(3,:,:)=trapmf(data.VEL,[-100,-99,vel.cloud(1),vel.cloud(2)]); % Houze 2014
m(4,:,:)=smf(data.TEMP,[temp.cloud(1),temp.cloud(2)]);

result(3,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(3,:,:),'CloudLiquid',plotIn);
end

if plotIn.plotMR
    plotMresult(data,m,result(4,:,:),'MixedPhase',plotIn);
end

%  Membership functions for large frozen
m=nan(4,size(data.DBZ,1),size(data.DBZ,2));
m(1,:,:)=smf(data.DBZ,[dbz.lfrozen(1),dbz.lfrozen(2)]); 
m(2,:,:)=trapmf(data.LDR,[ldr.lfrozen(1),ldr.lfrozen(2),ldr.lfrozen(3),ldr.lfrozen(4)]);
m(3,:,:)=smf(data.VEL,[vel.lfrozen(1),vel.lfrozen(2)]);
m(4,:,:)=trapmf(data.TEMP,[-100,-99,0,6]);

result(4,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(4,:,:),'LargeFrozen',plotIn);
end

%  Membership functions for small frozen
m=nan(4,size(data.DBZ,1),size(data.DBZ,2));
m(1,:,:)=trapmf(data.DBZ,[-100,-99,dbz.sfrozen(1),dbz.sfrozen(2)]);
m(2,:,:)=trapmf(data.LDR,[ldr.sfrozen(1),ldr.sfrozen(2),ldr.sfrozen(3),ldr.sfrozen(4)]);
m(3,:,:)=trapmf(data.VEL,[-100,-99,vel.sfrozen(1),vel.sfrozen(2)]);
m(4,:,:)=trapmf(data.TEMP,[-100,-99,temp.sfrozen(1),temp.sfrozen(2)]);

result(5,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(5,:,:),'SmallFrozen',plotIn);
end

clear m

[maxAll,classOut]=max(result,[],1,'omitnan');
maxAll=squeeze(maxAll);
classOut=squeeze(classOut);

% Remove nans
classOut(isnan(data.DBZ))=nan;

if plotIn.plotMax
    plotResMax(data,result,maxAll,plotIn);
end

classOut(classOut==5)=6;
classOut(classOut==4)=5;

classOut(isnan(dbzOrig))=nan;
end