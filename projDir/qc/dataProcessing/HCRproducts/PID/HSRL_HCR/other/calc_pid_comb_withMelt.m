function[classOut]=calc_pid_comb(data,plotIn)
%   Membership functions for particle detection

memCoeffsComb
% DBZ, HCR LDR, VEL, TEMP, BACKSCATTER, HSRL LDR
w=[22 16 16 16 14 16];
%w=[30 20 20 20 14 16];

% pid_hcr (number before post processing)
%  1 Rain
%  2 Drizzle
%  3 Cloud liquid
%  4 Mixed phase
%  5 Large frozen
%  6 Small frozen

result=nan(6,size(data.HCR_LDR,1),size(data.HCR_LDR,2));

%  Membership functions for rain
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
m(1,:,:)=smf(data.HCR_DBZ,[dbz.rain(1),dbz.rain(2)]);  % Rain
m(2,:,:)=zmf(data.HCR_LDR,[hcrldr.rain(1),hcrldr.rain(2)]);
m(3,:,:)=smf(data.HCR_VEL,[vel.rain(1),vel.rain(2)]); % Beard and Pruppacher 1969
m(4,:,:)=smf(data.TEMP,[temp.rain(1),temp.rain(2)]);
m(5,:,:)=smf(data.HSRL_Aerosol_Backscatter_Coefficient,[back.rain(1),back.rain(2)]);
m(6,:,:)=zmf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.rain(1),hsrlldr.rain(2)]);

result(1,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(1,:,:),'Rain',plotIn);
end

%  Membership functions for drizzle
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
m(1,:,:)=trapmf(data.HCR_DBZ,[dbz.drizzle(1),dbz.drizzle(2),dbz.drizzle(3),dbz.drizzle(4)]);
m(2,:,:)=zmf(data.HCR_LDR,[hcrldr.drizzle(1),hcrldr.drizzle(2)]);
m(3,:,:)=trapmf(data.HCR_VEL,[vel.drizzle(1),vel.drizzle(2),vel.drizzle(3),vel.drizzle(4)]);% Houze 2014, Beard and Pruppacher 1969
m(4,:,:)=smf(data.TEMP,[temp.drizzle(1),temp.drizzle(2)]);
m(5,:,:)=smf(data.HSRL_Aerosol_Backscatter_Coefficient,[back.drizzle(1),back.drizzle(2)]);
m(6,:,:)=zmf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.drizzle(1),hsrlldr.drizzle(2)]);

result(2,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(2,:,:),'Drizzle',plotIn);
end

%  Membership functions for cloud liquid
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
m(1,:,:)=zmf(data.HCR_DBZ,[dbz.cloud(1),dbz.cloud(2)]);
m(2,:,:)=zmf(data.HCR_LDR,[hcrldr.cloud(1),hcrldr.cloud(2)]);
m(3,:,:)=zmf(data.HCR_VEL,[vel.cloud(1),vel.cloud(2)]); % Houze 2014
m(4,:,:)=smf(data.TEMP,[temp.cloud(1),temp.cloud(2)]);
m(5,:,:)=smf(data.HSRL_Aerosol_Backscatter_Coefficient,[back.cloud(1),back.cloud(2)]);
m(6,:,:)=zmf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.cloud(1),hsrlldr.cloud(2)]);

result(3,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(3,:,:),'CloudLiquid',plotIn);
end

%  Membership functions for mixed phase
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
%m(1,:,:)=smf(data.HCR_DBZ,[dbz.mixed(1),dbz.mixed(2)]);
m(1,:,:)=1;
m(2,:,:)=trapmf(data.HCR_LDR,[hcrldr.mixed(1),hcrldr.mixed(2),hcrldr.mixed(3),hcrldr.mixed(4)]);
m(3,:,:)=smf(data.HCR_VEL,[vel.mixed(1),vel.mixed(2)]);
m(4,:,:)=trapmf(data.TEMP,[temp.mixed(1),temp.mixed(2),temp.mixed(3),temp.mixed(4)]);
m(5,:,:)=0;
m(6,:,:)=zmf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.mixed(1),hsrlldr.mixed(2)]);

result(4,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(4,:,:),'MixedPhase',plotIn);
end

%  Membership functions for large frozen
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
m(1,:,:)=smf(data.HCR_DBZ,[dbz.lfrozen(1),dbz.lfrozen(2)]); 
m(2,:,:)=trapmf(data.HCR_LDR,[hcrldr.lfrozen(1),hcrldr.lfrozen(2),hcrldr.lfrozen(3),hcrldr.lfrozen(4)]);
m(3,:,:)=smf(data.HCR_VEL,[vel.lfrozen(1),vel.lfrozen(2)]);
m(4,:,:)=zmf(data.TEMP,[0,6]);
m(5,:,:)=zmf(data.HSRL_Aerosol_Backscatter_Coefficient,[back.lfrozen(1),back.lfrozen(2)]);
m(6,:,:)=smf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.lfrozen(1),hsrlldr.lfrozen(2)]);

result(5,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(5,:,:),'LargeFrozen',plotIn);
end

%  Membership functions for small frozen
m=nan(6,size(data.HCR_DBZ,1),size(data.HCR_DBZ,2));
m(1,:,:)=zmf(data.HCR_DBZ,[dbz.sfrozen(1),dbz.sfrozen(2)]);
m(2,:,:)=trapmf(data.HCR_LDR,[hcrldr.sfrozen(1),hcrldr.sfrozen(2),hcrldr.sfrozen(3),hcrldr.sfrozen(4)]);
m(3,:,:)=zmf(data.HCR_VEL,[vel.sfrozen(1),vel.sfrozen(2)]);
m(4,:,:)=zmf(data.TEMP,[temp.sfrozen(1),temp.sfrozen(2)]);
m(5,:,:)=zmf(data.HSRL_Aerosol_Backscatter_Coefficient,[back.sfrozen(1),back.sfrozen(2)]);
m(6,:,:)=smf(data.HSRL_Particle_Linear_Depolarization_Ratio,[hsrlldr.sfrozen(1),hsrlldr.sfrozen(2)]);

result(6,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(6,:,:),'SmallFrozen',plotIn);
end

clear m

[maxAll,classOut]=max(result,[],1,'omitnan');
maxAll=squeeze(maxAll);
classOut=squeeze(classOut);

% Remove nans
classOut(isnan(data.HCR_DBZ))=nan;

if plotIn.plotMax
    plotResMax(data,result,maxAll,plotIn);
end

end