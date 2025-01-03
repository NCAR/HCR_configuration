function[classOut]=calc_pid_hcr_postProcess(dBZ,data,postProcess)

%data.HCR_VEL=abs(data.HCR_VEL);
data.HCR_VEL(data.elevation>0)=-data.HCR_VEL(data.elevation>0);

%   Membership functions for particle detection
% 1:Beta  2:Delta

w=[40 15 10 10 25];%w=[30 15 15 20 20];

% pid_hcr
%  1 no signal -> nan after final reassignment
%  2 cloud liquid -> 1 after final reassignment
%  3 Drizzle -> 2 after final reassignment
%  4 Rain -> 3 after final reassignment
%  5 SLW -> 4 after final reassignment
%  6 Ice crystals -> 5 after final reassignment
%  7 Snow -> 6 after final reassignment
%  8 wet snow/rimed ice -> 7 after final reassignment

result=nan(8,size(data.HCR_LDR,1),size(data.HCR_LDR,2));

%  Membership functions for no signal
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)= zmf(dBZ,[-30,-32]);
m(2,:,:)=zmf(data.HCR_LDR,[-25, -27]);
m(3,:,:)=trapmf(data.HCR_VEL,[-6,-4,4,6]);
m(4,:,:)=zmf(data.HCR_WIDTH,[0.2, 0.1]);
m(5,:,:)=smf(data.temp,[234., 300.]);

result(1,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for cloud liquid
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)= trapmf(dBZ,[-32,-30,-17,-15]); % cloud liquid
m(2,:,:)=zmf(data.HCR_LDR,[-25, -27]);
m(3,:,:)=trapmf(data.HCR_VEL,[-6,-5,1.0, 2.0]);
%m(3,:,:)=trapmf(data.HCR_VEL,[-6,-5,-0.5,-0.2]);
%m(3,:,:)=trapmf(data.HCR_VEL,[-1,0,0.5,1]);
m(4,:,:)=zmf(data.HCR_WIDTH,[0.2, 0.1]);
%m(4,:,:)=zmf(data.HCR_WIDTH,[0.4, 0.3]);
m(5,:,:)=smf(data.temp,[272., 275.]);

result(2,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for drizzle
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-17,-14,5,8]);
m(2,:,:)=zmf(data.HCR_LDR,[-25, -27]);
m(3,:,:)=trapmf(data.HCR_VEL,[0,0.5,1,2]);
m(4,:,:)=zmf(data.HCR_WIDTH,[0.2, 0.1]);
%m(4,:,:)=zmf(data.HCR_WIDTH,[0.4, 0.3]);
m(5,:,:)=smf(data.temp,[246., 248.]);

result(3,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for rain
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=smf(dBZ,[3,5]);  % Rain
m(2,:,:)= zmf(data.HCR_LDR,[-22, -27]);
m(3,:,:)=trapmf(data.HCR_VEL,[2,3,6,8]);
m(4,:,:)= smf(data.HCR_WIDTH,[0.2, 0.3]);
%m(4,:,:)= smf(data.HCR_WIDTH,[0.4, 0.5]);
m(5,:,:)=smf(data.temp,[272., 275.]);

result(4,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for SLW
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=zmf(dBZ,[-14,-17]);
m(2,:,:)= zmf(data.HCR_LDR,[-25, -27]);
m(3,:,:)=trapmf(data.HCR_VEL,[-6,-5,1.0, 2.0]);
%m(3,:,:)=trapmf(data.HCR_VEL,[-6,-5,-0.5,-0.2]);
%m(3,:,:)=trapmf(data.HCR_VEL,[-1,0,1,2]);
m(4,:,:)=zmf(data.HCR_WIDTH,[0.2, 0.1]);
%m(4,:,:)=zmf(data.HCR_WIDTH,[0.4, 0.3]);
m(5,:,:)=zmf(data.temp,[273,271]);

result(5,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for ice
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-25,-20,5,10]);
%m(1,:,:)=trapmf(dBZ,[-32,-30,5,10]);
%m(1,:,:)=trapmf(dBZ,[-25,-20,5,10]);
m(2,:,:)=trapmf(data.HCR_LDR,[-25, -22,-15,-12]);
%m(2,:,:)=trapmf(data.HCR_LDR,[-25, -22,-20,-18]);
m(3,:,:)=trapmf(data.HCR_VEL,[-1,0,1,4]);
%m(3,:,:)=trapmf(data.HCR_VEL,[-1,0,1,2]);
m(4,:,:)=zmf(data.HCR_WIDTH,[0.7, 0.9]);
%m(4,:,:)=zmf(data.HCR_WIDTH,[0.4, 0.3]);
m(5,:,:)=zmf(data.temp,[273,271]);

result(6,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for snow
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[10,12,18,20]); % Snow
m(2,:,:)=trapmf(data.HCR_LDR,[-22,-18,-16, -14]);
m(3,:,:)=trapmf(data.HCR_VEL,[0.8,1.0,1.2,1.4]);
m(4,:,:)= smf(data.HCR_WIDTH,[0.2, 0.3]);
%m(4,:,:)= smf(data.HCR_WIDTH,[0.4, 0.5]);
m(5,:,:)=zmf(data.temp,[273,271]);

result(7,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

%  Membership functions for snow/rimed ice
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-3,-1,10,15]);
m(2,:,:)=trapmf(data.HCR_LDR,[-20, -17,-8,-6]);
m(3,:,:)=trapmf(data.HCR_VEL,[0.5,1.0,3,4]);
m(4,:,:)= smf(data.HCR_WIDTH,[0.2, 0.3]);
%m(4,:,:)=trapmf(data.HCR_WIDTH,[0.3,0.4,0.6,0.7]);
m(5,:,:)=trapmf(data.temp,[268,270,275,277]);

result(8,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
        +m(4,:,:)*w(4)+m(5,:,:)*w(5);

clear temp sigmav vel dBZ

clear m

inoldr=find(isnan(data.HCR_LDR)==1);
result(8,inoldr)=0;

max1=squeeze(nanmax(result,[],1));
maxMat=repmat(max1,[1,1,8]);
maxMat=permute(maxMat,[3,1,2]);
resMinusMax=result-maxMat;
zerosMat=zeros(size(result));
zerosMat(resMinusMax==0)=1;

classOut=nan(size(data.HCR_LDR));

for ii=1:8
    testMat=squeeze(zerosMat(ii,:,:));
    classOut(isnan(classOut) & testMat==1)=ii;
end

clearvars -except data classOut postProcess
%% Post processing
if postProcess
    meltLayer=data.MELTING_LAYER;
    meltLayer(~isnan(meltLayer) & meltLayer<20)=10;
    meltLayer(~isnan(meltLayer) & meltLayer>=20)=20;
    
    % No frozen particles in strong downward motion
    classOut(meltLayer==10 & data.HCR_VEL>2.5 & ...
        (classOut==6 | classOut==7 | classOut==8))=4;
    
    % No frozen precipitation in warm region
    replaceMat=zeros(size(classOut));
    replaceMat(meltLayer==10 & (classOut==6 | classOut==7))=1;
    
    % Replace with closest warm pixel
    [oldR oldC]=find(~isnan(classOut) & replaceMat==0 & meltLayer==10);
    [addR addC]=find(replaceMat==1);
    idx = knnsearch([oldR oldC], [addR addC]);
    nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
    classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;
    
    % Updrafts have no rain, no drizzle, and no snow.
    replaceMat=zeros(size(classOut));
    replaceMat(data.HCR_VEL<0 & (classOut==3 | classOut==4 | classOut==7))=1;
    
    % Replace with closest pixel
    [oldR oldC]=find(~isnan(classOut) & replaceMat==0);
    [addR addC]=find(replaceMat==1);
    idx = knnsearch([oldR oldC], [addR addC]);
    nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
    classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;
end
%% Reassign values
classOut(classOut==1)=nan;
classOut=classOut-1;
end