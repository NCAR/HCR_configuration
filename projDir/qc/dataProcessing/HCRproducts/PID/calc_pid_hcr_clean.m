function[classOut, m]=calc_pid_hcr_clean(dBZ,LDR, vel,sigmav,temp)
vel=abs(vel);
%   Membership functions for particle detection
% 1:Beta  2:Delta

w=[40 15 15 10 20];

% pid_hcr
%  1 no signal
%  2 cloud liquid
%  3 Drizzle
%  4 Rain
%  5 SLW
%  6 Ice crystals
%  7 Snow
%  8 wet snow/rimed ice

%  Membership functions for Reflectivity
m=nan(8,5,size(dBZ,1),size(dBZ,2));
m(1,1,:,:)= zmf(dBZ,[-30,-32]); % no signal
m(2,1,:,:)= trapmf(dBZ,[-32,-30,-17,-15]); % cloud liquid
m(3,1,:,:)=trapmf(dBZ,[-17,-14,5,8]);  % Drizzle
m(4,1,:,:)=smf(dBZ,[3,5]);  % Rain
m(5,1,:,:)=zmf(dBZ,[-14,-17]);   % SLW
m(6,1,:,:)=trapmf(dBZ,[-25,-20,5,10]);  % Ice
m(7,1,:,:)=trapmf(dBZ,[10,12,18,20]); % Snow
m(8,1,:,:)=trapmf(dBZ,[-3,-1,10,15]); % wet snow/rimed ice

%  Membership functions for LDR
m(1,2,:,:)=zmf(LDR,[-25, -27]);% no signal
m(2,2,:,:)=zmf(LDR,[-25, -27]); % cloud liquid
m(3,2,:,:)=zmf(LDR,[-25, -27]); % Drizzle
m(4,2,:,:)= zmf(LDR,[-22, -27]); % Rain
m(5,2,:,:)= zmf(LDR,[-25, -27]);  % SLW
m(6,2,:,:)=trapmf(LDR,[-25, -22,-20,-18]); % Ice
m(7,2,:,:)=trapmf(LDR,[-22,-18,-16, -14]); % Snow
m(8,2,:,:)=trapmf(LDR,[-20, -17,-8,-6]); % wet snow/rimed ice

%  Membership functions for velocity
m(1,3,:,:)=trapmf(vel,[-6,-4,4,6]);% no signal
m(2,3,:,:)=trapmf(vel,[-1,0,0.5,1]); % cloud liquid
m(3,3,:,:)=trapmf(vel,[0,0.5,1,2]); % Drizzle
m(4,3,:,:)=trapmf(vel,[2,3,6,8]); % Rain
m(5,3,:,:)=trapmf(vel,[-1,0,1,2]);  % SLW
m(6,3,:,:)=trapmf(vel,[-1,0,1,2]); % Ice
m(7,3,:,:)=trapmf(vel,[0.8,1.0,1.2,1.4]); % Snow
m(8,3,:,:)=trapmf(vel,[0.5,1.0,3,4]); % wet snow/rimed ice

%  Membership functions for sigmav
m(1,4,:,:)=zmf(sigmav,[0.2, 0.1]);% no signal
m(2,4,:,:)=zmf(sigmav,[0.4, 0.3]); % cloud liquid
m(3,4,:,:)=zmf(sigmav,[0.4, 0.3]); % Drizzle
m(4,4,:,:)= smf(sigmav,[0.4, 0.5]); % Rain
m(5,4,:,:)= zmf(sigmav,[0.4, 0.3]);  % SLW
m(6,4,:,:)=zmf(sigmav,[0.4, 0.3]); % Ice
m(7,4,:,:)=smf(sigmav,[0.4, 0.5]); % Snow
m(8,4,:,:)=trapmf(sigmav,[0.3,0.4,0.6,0.7]); % wet snow/rimed ice

%  Membership functions for Temperature
m(1,5,:,:)=smf(temp,[234., 300.]);% no signal
m(2,5,:,:)=smf(temp,[272., 275.]); % cloud liquid
m(3,5,:,:)=smf(temp,[246., 248.]); % Drizzle
m(4,5,:,:)=smf(temp,[272., 275.]); % Rain
m(5,5,:,:)=zmf(temp,[273,271]);  % SLW
m(6,5,:,:)=zmf(temp,[273,271]); % Ice
m(7,5,:,:)=zmf(temp,[273,271]); % Snow
m(8,5,:,:)=trapmf(temp,[268,270,275,277]); % wet snow/rimed ice

%  Sum the weights
result=nan(8,size(dBZ,1),size(dBZ,2));

for ii=1:8
    mii=squeeze(m(ii,:,:,:));
    result(ii,:,:)=mii(1,:,:)*w(1)+mii(2,:,:)*w(2)+mii(3,:,:)*w(3)...
        +mii(4,:,:)*w(4)+mii(5,:,:)*w(5);
end

inoldr=find(isnan(LDR)==1);
result(8,inoldr)=0;

max1=squeeze(nanmax(result,[],1));
maxMat=repmat(max1,[1,1,8]);
maxMat=permute(maxMat,[3,1,2]);
resMinusMax=result-maxMat;
zerosMat=zeros(size(result));
zerosMat(resMinusMax==0)=1;

classOut=nan(size(dBZ));

for ii=1:8
    testMat=squeeze(zerosMat(ii,:,:));
    classOut(isnan(classOut) & testMat==1)=ii;
end

end