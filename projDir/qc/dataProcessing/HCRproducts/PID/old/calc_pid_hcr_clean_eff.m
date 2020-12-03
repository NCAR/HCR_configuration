function[classOut]=calc_pid_hcr_clean_eff(dBZ,LDR, vel,sigmav,temp)
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

result=nan(8,size(LDR,1),size(LDR,2));

%  Membership functions for no signal
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)= zmf(dBZ,[-30,-32]); % no signal
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[-6,-4,4,6]);
m(4,:,:)=zmf(sigmav,[0.2, 0.1]);
m(5,:,:)=smf(temp,[234., 300.]);

result(1,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for cloud liquid
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)= trapmf(dBZ,[-32,-30,-17,-15]); % cloud liquid
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[-1,0,0.5,1]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=smf(temp,[272., 275.]);

result(2,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for drizzle
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-17,-14,5,8]);  % Drizzle
m(2,:,:)=zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[0,0.5,1,2]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=smf(temp,[246., 248.]);

result(3,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for rain
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=smf(dBZ,[3,5]);  % Rain
m(2,:,:)= zmf(LDR,[-22, -27]);
m(3,:,:)=trapmf(vel,[2,3,6,8]);
m(4,:,:)= smf(sigmav,[0.4, 0.5]);
m(5,:,:)=smf(temp,[272., 275.]);

result(4,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for SLW
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=zmf(dBZ,[-14,-17]);   % SLW
m(2,:,:)= zmf(LDR,[-25, -27]);
m(3,:,:)=trapmf(vel,[-1,0,1,2]);
m(4,:,:)= zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=zmf(temp,[273,271]);

result(5,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for ice
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-25,-20,5,10]);  % Ice
m(2,:,:)=trapmf(LDR,[-25, -22,-20,-18]);
m(3,:,:)=trapmf(vel,[-1,0,1,2]);
m(4,:,:)=zmf(sigmav,[0.4, 0.3]);
m(5,:,:)=zmf(temp,[273,271]);

result(6,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for snow
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[10,12,18,20]); % Snow
m(2,:,:)=trapmf(LDR,[-22,-18,-16, -14]);
m(3,:,:)=trapmf(vel,[0.8,1.0,1.2,1.4]);
m(4,:,:)=smf(sigmav,[0.4, 0.5]);
m(5,:,:)=zmf(temp,[273,271]);

result(7,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for snow/rimed ice
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-3,-1,10,15]); % wet snow/rimed ice
m(2,:,:)=trapmf(LDR,[-20, -17,-8,-6]); % wet snow/rimed ice
m(3,:,:)=trapmf(vel,[0.5,1.0,3,4]); % wet snow/rimed ice
m(4,:,:)=trapmf(sigmav,[0.3,0.4,0.6,0.7]); % wet snow/rimed ice
m(5,:,:)=trapmf(temp,[268,270,275,277]); % wet snow/rimed ice

result(8,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

clear temp sigmav vel dBZ

clear m

inoldr=find(isnan(LDR)==1);
result(8,inoldr)=0;

max1=squeeze(nanmax(result,[],1));
maxMat=repmat(max1,[1,1,8]);
maxMat=permute(maxMat,[3,1,2]);
resMinusMax=result-maxMat;
zerosMat=zeros(size(result));
zerosMat(resMinusMax==0)=1;

classOut=nan(size(LDR));

for ii=1:8
    testMat=squeeze(zerosMat(ii,:,:));
    classOut(isnan(classOut) & testMat==1)=ii;
end

end