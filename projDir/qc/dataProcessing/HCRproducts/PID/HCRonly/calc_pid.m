function[classOut]=calc_pid(dBZ,data,postProcess,plotIn)

data.VEL_MASKED(data.elevation>0)=-data.VEL_MASKED(data.elevation>0);

%   Membership functions for particle detection
% 1:Beta  2:Delta

w=[40 15 10 10 25];%w=[30 15 15 20 20];

% pid_hcr (number before post processing)
%  1 Rain (1)
%  2 Supercooled rain (post-processing)
%  3 Drizzle (2)
%  4 Supercooled drizzle (post-processing)
%  5 Cloud liquid (3)
%  6 Supercooled cloud liquid (post processing)
%  7 Mixed phase (4)
%  8 Large frozen (5)
%  9 Small frozen (6)

result=nan(6,size(data.LDR,1),size(data.LDR,2));

%  Membership functions for rain
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=smf(dBZ,[3,5]);  % Rain
m(2,:,:)=zmf(data.LDR,[-27,-22]);
m(3,:,:)=smf(data.VEL_MASKED,[2,3]);
m(4,:,:)=smf(data.WIDTH,[0.1,0.2]);
m(5,:,:)=smf(data.TEMP,[-2,2]);

result(1,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

if plotIn.plotMR
    plotMresult(data,m,result(1,:,:),'Rain',plotIn);
end

%  Membership functions for drizzle
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-18,-16,5,8]);
m(2,:,:)=zmf(data.LDR,[-27,-25]);
m(3,:,:)=trapmf(data.VEL_MASKED,[0,0.5,1,2]);
m(4,:,:)=zmf(data.WIDTH,[0.2,0.3]);
m(5,:,:)=smf(data.TEMP,[-50,-39]);

result(2,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

if plotIn.plotMR
    plotMresult(data,m,result(2,:,:),'Drizzle',plotIn);
end

%  Membership functions for cloud liquid
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=zmf(dBZ,[-16,-14]);
m(2,:,:)=zmf(data.LDR,[-27,-25]);
m(3,:,:)=zmf(data.VEL_MASKED,[1,2]);
m(4,:,:)=zmf(data.WIDTH,[0.1,0.2]);
m(5,:,:)=1;

result(3,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

if plotIn.plotMR
    plotMresult(data,m,result(3,:,:),'CloudLiquid',plotIn);
end

%  Membership functions for mixed phase
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-3,-1,20,25]);
m(2,:,:)=trapmf(data.LDR,[-20,-17,-8,-6]);
m(3,:,:)=trapmf(data.VEL_MASKED,[0.5,1,3,4]);
m(4,:,:)=smf(data.WIDTH,[0.2, 0.3]);
m(5,:,:)=trapmf(data.TEMP,[-2,0,3,6]);

result(4,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

inoldr=find(isnan(data.LDR));
result(4,inoldr)=0;

if plotIn.plotMR
    plotMresult(data,m,result(4,:,:),'MixedPhase',plotIn);
end

%  Membership functions for large frozen
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[7,9,18,20]); 
m(2,:,:)=trapmf(data.LDR,[-22,-20,-16,-14]);
m(3,:,:)=trapmf(data.VEL_MASKED,[0.8,1,2.5,3.5]);
m(4,:,:)=smf(data.WIDTH,[0.2, 0.3]);
m(5,:,:)=zmf(data.TEMP,[0,6]);

result(5,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

if plotIn.plotMR
    plotMresult(data,m,result(5,:,:),'LargeFrozen',plotIn);
end

%  Membership functions for small frozen
m=nan(5,size(dBZ,1),size(dBZ,2));
m(1,:,:)=trapmf(dBZ,[-25,-20,9,11]);
m(2,:,:)=trapmf(data.LDR,[-26, -24,-15,-12]);
m(3,:,:)=trapmf(data.VEL_MASKED,[-1,0,1,2]);
m(4,:,:)=zmf(data.WIDTH,[0.7,0.9]);
m(5,:,:)=zmf(data.TEMP,[-1,5]);

result(6,:,:)=m(1,:,:)*w(1)+m(2,:,:)*w(2)+m(3,:,:)*w(3)...
    +m(4,:,:)*w(4)+m(5,:,:)*w(5);

if plotIn.plotMR
    plotMresult(data,m,result(6,:,:),'SmallFrozen',plotIn);
end

clear m

max1=squeeze(nanmax(result,[],1));
maxMat=repmat(max1,[1,1,6]);
maxMat=permute(maxMat,[3,1,2]);
resMinusMax=result-maxMat;
zerosMat=zeros(size(result));
zerosMat(resMinusMax==0)=1;

classOut=nan(size(data.LDR));

for ii=1:6
    testMat=squeeze(zerosMat(ii,:,:));
    classOut(isnan(classOut) & testMat==1)=ii;
end

if plotIn.plotMax
    plotResMax(data,result,max1,plotIn);
end

clearvars -except data classOut postProcess

%% Add supercooled
classOut(classOut==6)=9;
classOut(classOut==5)=8;
classOut(classOut==4)=7;
classOut(classOut==3)=5;
classOut(classOut==2)=3;

% Supercooled rain
classOut(classOut==1 & data.TEMP<0)=2;

% Supercooled drizzle
classOut(classOut==3 & data.TEMP<0)=4;

% Supercooled cloud liquid
classOut(classOut==5 & data.TEMP<0)=6;

%% Post processing
if postProcess
    meltLayer=data.MELTING_LAYER;
    meltLayer(~isnan(meltLayer) & meltLayer<20)=10;
    meltLayer(~isnan(meltLayer) & meltLayer>=20)=20;
    
    % No small frozen in warm region
    replaceMat=zeros(size(classOut));
    replaceMat(meltLayer==10 & classOut==9)=1;
        
    % Replace with closest warm pixel
    [oldR oldC]=find(~isnan(classOut) & replaceMat==0 & meltLayer==10);
    [addR addC]=find(replaceMat==1);
    idx = knnsearch([oldR oldC], [addR addC]);
    nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
    classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;
    
end
%% Remove nans
classOut(isnan(data.DBZ_MASKED))=nan;
end