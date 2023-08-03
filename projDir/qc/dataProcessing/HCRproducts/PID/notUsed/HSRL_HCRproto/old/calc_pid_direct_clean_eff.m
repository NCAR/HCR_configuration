function[classOut, m]=calc_pid_direct_clean(Beta, Delta,dBZ,LDR, vel, sigmav,temp)
vel=abs(vel);
%   Membership functions for particle detection

% 1:Beta  2:Delta
w=[25 15 25 10 5 5 15];

%   w=[20 15 25 10 5 5 20];

% pid_Precip_Aerosol
%  1 no signal
%  2 cloud liquid
%  3 Drizzle
%  4 Rain
%  5 SLW
%  6 Ice crystals
%  7 Snow
%  8 wet snow/rimed ice
%  9 Aerosol

result=nan(9,size(dBZ,1),size(dBZ,2));

%  Membership functions for no signal
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)= zmf(dBZ,[-30,-32]); % no signal
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=zmf(vel,[2 0]);
m(4,:,:)=zmf(sigmav,[0.2, 0.1]);
m(5,:,:)=smf(temp,[234., 236.]);
m(6,:,:)= zmf(Beta,[0.4e-10,0.2e-9]);
m(7,:,:)=smf(Delta,[0, 0.01]);

result(1,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for cloud liquid
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)= trapmf(dBZ,[-32,-30,-17,-15]); % cloud liquid
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[-1,0,0.5,1]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=smf(temp,[273., 275.]);
m(6,:,:)=trapmf(Beta,[1.8e-5,2e-5,8e-4,1.0e-3]);
m(7,:,:)=trapmf(Delta,[-0.1, 0.01,0.04,0.05]);

result(2,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for drizzle
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-10,-8,5,8]);  % Drizzle
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[0,0.5,1,2]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=smf(temp,[246., 248.]);
m(6,:,:)=trapmf(Beta,[0.8e-6,1.2e-6,0.8e-4,1.0e-4]);
m(7,:,:)=trapmf(Delta,[0.04,0.05,0.23, 0.25]);

result(3,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for rain
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=smf(dBZ,[3,5]);  % Rain
m(2,:,:)= zmf(LDR,[-22, -27]);
m(3,:,:)=trapmf(vel,[2,3,6,8]);
m(4,:,:)= smf(sigmav,[0.4, 0.5]);
m(5,:,:)=smf(temp,[272., 275.]);
m(6,:,:)=trapmf(Beta,[0.8e-6,1.2e-6,0.8e-4,1.0e-4]);
m(7,:,:)=trapmf(Delta,[0.04,0.05,0.23, 0.25]);

result(4,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for SLW
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=zmf(dBZ,[-14,-17]);   % SLW
m(2,:,:)= zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[-1,0,0.5,1]);
m(4,:,:)= zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=zmf(temp,[273,271]);
m(6,:,:)=trapmf(Beta,[1.8e-5,2e-5,8e-4,1.0e-3]);
m(7,:,:)=trapmf(Delta,[-0.1, 0.01,0.04,0.05]);

result(5,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for ice crystals
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-25,-20,5,10]);  % Ice crystal
m(2,:,:)=trapmf(LDR,[-25, -22,-20,-18]);
m(3,:,:)=trapmf(vel,[-1,0,1,2]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=zmf(temp,[273,271]);
m(6,:,:)=trapmf(Beta,[0.8e-6,1.2e-6,0.8e-4,1.0e-4]);
m(7,:,:)= trapmf(Delta,[0.23,0.25,0.50, 0.6]);

result(6,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for snow
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[5,10,15,20]); % Snow
m(2,:,:)=trapmf(LDR,[-22,-18,-16, -14]);
m(3,:,:)=trapmf(vel,[0,0.5,1,2]);
m(4,:,:)=smf(sigmav,[0.4, 0.5]);
m(5,:,:)=zmf(temp,[273,271]);
m(6,:,:)=trapmf(Beta,[0.8e-6,1.2e-6,0.8e-4,1.0e-4]);
m(7,:,:)= trapmf(Delta,[0.23,0.25,0.50, 0.6]);

result(7,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for wet snow/rimed ice
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-3,-1,10,15]); % wet snow/rimed ice
m(2,:,:)=trapmf(LDR,[-20, -17,-8,-6]);
m(3,:,:)=trapmf(vel,[0.5,1.0,3,4]);
m(4,:,:)=trapmf(sigmav,[0.3,0.4,0.6,0.7]);
m(5,:,:)=trapmf(temp,[248,250,275,277]);
m(6,:,:)=trapmf(Beta,[0.8e-6,1.2e-6,0.8e-4,1.0e-4]);
m(7,:,:)= trapmf(Delta,[0.23,0.25,0.50, 0.6]);

result(8,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

%  Membership functions for aerosol
m=nan(7,size(dBZ,1),size(dBZ,2));
m(1,:,:)= zmf(dBZ,[-30,-32]); % Aerosol; this is same as no signal
m(2,:,:)=zmf(LDR,[-25, -27]);% Aerosol
m(3,:,:)=zmf(vel,[2 0]);% Aerosol
m(4,:,:)=zmf(sigmav,[0.2, 0.1]);% Aeerosol
m(5,:,:)=smf(temp,[234., 236.]);% Aerosol
m(6,:,:)=trapmf(Beta,[0.2e-9,0.4e-9,1.8e-5,2e-5]);  % Aerosol
m(7,:,:)=smf(Delta,[0, 0.01]); % Aerosol

result(9,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5)+m(6,:,:)*w(6)+m(7,:,:)*w(7);

clear m

inoldr=find(isnan(LDR));
result(8,inoldr)=0;

max1=squeeze(nanmax(result,[],1));
maxMat=repmat(max1,[1,1,9]);
maxMat=permute(maxMat,[3,1,2]);
resMinusMax=result-maxMat;
zerosMat=zeros(size(result));
zerosMat(resMinusMax==0)=1;

classOut=nan(size(dBZ));

for ii=1:9
    testMat=squeeze(zerosMat(ii,:,:));
    classOut(isnan(classOut) & testMat==1)=ii;
end

end