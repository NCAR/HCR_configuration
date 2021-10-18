function[classOut]=calc_pid_ldr(DBZ,data,postProcess,plotIn)

% Make all down motion positive
data.VEL_MASKED(data.elevation>0)=-data.VEL_MASKED(data.elevation>0);

% Remove data with too much convectivity
data.VEL_MASKED(data.CONVECTIVITY>0.4)=nan;

%   Membership functions for particle detection
% 1:Beta  2:Delta

memCoeffsLDR
% DBZ, LDR, VEL, WIDTH, TEMP
w=[30 20 20 0 30];%w=[30 15 15 20 20];

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
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=smf(DBZ,[dbz.rain(1),dbz.rain(2)]);  % Rain
m(2,:,:)=zmf(data.LDR,[ldr.rain(1),ldr.rain(2)]);
m(3,:,:)=smf(data.VEL_MASKED,[vel.rain(1),vel.rain(2)]); % Beard and Pruppacher 1969
m(4,:,:)=smf(data.WIDTH,[width.rain(1),width.rain(2)]);
m(5,:,:)=1;

result(1,:,:)=sum(m.*w',1);

% Adjust weights
% mNan=nan(size(m));
% mTest=m(1,:,:);
% mTest(isnan(DBZ))=nan;
% mNan(1,:,:)=mTest;
% mTest=m(2,:,:);
% mTest(isnan(data.LDR))=nan;
% mNan(2,:,:)=mTest;
% mTest=m(3,:,:);
% mTest(isnan(data.VEL_MASKED))=nan;
% mNan(3,:,:)=mTest;
% mTest=m(4,:,:);
% mTest(isnan(data.WIDTH))=nan;
% mNan(4,:,:)=mTest;
% mTest=m(5,:,:);
% mTest(isnan(data.TEMP))=nan;
% mNan(5,:,:)=mTest;
% 
% W = bsxfun(@times,~isnan(mNan),w');
% result1=squeeze(sum(m.*W,1)./sum(W,1));

if plotIn.plotMR
    plotMresult(data,m,result(1,:,:),'Rain',plotIn);
end

%  Membership functions for drizzle
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=trapmf(DBZ,[dbz.drizzle(1),dbz.drizzle(2),dbz.drizzle(3),dbz.drizzle(4)]);
m(2,:,:)=zmf(data.LDR,[ldr.drizzle(1),ldr.drizzle(2)]);
m(3,:,:)=trapmf(data.VEL_MASKED,[vel.drizzle(1),vel.drizzle(2),vel.drizzle(3),vel.drizzle(4)]);% Houze 2014, Beard and Pruppacher 1969
m(4,:,:)=zmf(data.WIDTH,[width.drizzle(1),width.drizzle(2)]);
m(5,:,:)=1;

result(2,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(2,:,:),'Drizzle',plotIn);
end

%  Membership functions for cloud liquid
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=zmf(DBZ,[dbz.cloud(1),dbz.cloud(2)]);
m(2,:,:)=zmf(data.LDR,[ldr.cloud(1),ldr.cloud(2)]);
m(3,:,:)=zmf(data.VEL_MASKED,[vel.cloud(1),vel.cloud(2)]); % Houze 2014
m(4,:,:)=zmf(data.WIDTH,[width.cloud(1),width.cloud(2)]);
m(5,:,:)=1;

result(3,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(3,:,:),'CloudLiquid',plotIn);
end

%  Membership functions for mixed phase
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=smf(DBZ,[dbz.mixed(1),dbz.mixed(2)]);
m(2,:,:)=trapmf(data.LDR,[ldr.mixed(1),ldr.mixed(2),ldr.mixed(3),ldr.mixed(4)]);
m(3,:,:)=smf(data.VEL_MASKED,[vel.mixed(1),vel.mixed(2)]);
m(4,:,:)=smf(data.WIDTH,[width.mixed(1),width.mixed(2)]);
m(5,:,:)=trapmf(data.TEMP,[temp.mixed(1),temp.mixed(2),temp.mixed(3),temp.mixed(4)]);

result(4,:,:)=sum(m.*w',1);

% inoldr=find(isnan(data.LDR));
% result(4,inoldr)=0;

if plotIn.plotMR
    plotMresult(data,m,result(4,:,:),'MixedPhase',plotIn);
end

%  Membership functions for large frozen
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=smf(DBZ,[dbz.lfrozen(1),dbz.lfrozen(2)]); 
m(2,:,:)=trapmf(data.LDR,[ldr.lfrozen(1),ldr.lfrozen(2),ldr.lfrozen(3),ldr.lfrozen(4)]);
m(3,:,:)=smf(data.VEL_MASKED,[vel.lfrozen(1),vel.lfrozen(2)]);
m(4,:,:)=smf(data.WIDTH,[width.lfrozen(1),width.lfrozen(2)]);
m(5,:,:)=zmf(data.TEMP,[0,6]);

result(5,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(5,:,:),'LargeFrozen',plotIn);
end

%  Membership functions for small frozen
m=nan(5,size(DBZ,1),size(DBZ,2));
m(1,:,:)=zmf(DBZ,[dbz.sfrozen(1),dbz.sfrozen(2)]);
m(2,:,:)=trapmf(data.LDR,[ldr.sfrozen(1),ldr.sfrozen(2),ldr.sfrozen(3),ldr.sfrozen(4)]);
m(3,:,:)=trapmf(data.VEL_MASKED,[vel.sfrozen(1),vel.sfrozen(2),vel.sfrozen(3),vel.sfrozen(4)]);
m(4,:,:)=zmf(data.WIDTH,[width.sfrozen(1),width.sfrozen(2)]);
m(5,:,:)=zmf(data.TEMP,[temp.sfrozen(1),temp.sfrozen(2)]);

result(6,:,:)=sum(m.*w',1);

if plotIn.plotMR
    plotMresult(data,m,result(6,:,:),'SmallFrozen',plotIn);
end

clear m

[maxAll,classOut]=max(result,[],1,'omitnan');
maxAll=squeeze(maxAll);
classOut=squeeze(classOut);

% Remove nans
classOut(isnan(data.DBZ_MASKED))=nan;

% max1=squeeze(nanmax(result,[],1));
% maxMat=repmat(max1,[1,1,6]);
% maxMat=permute(maxMat,[3,1,2]);
% resMinusMax=result-maxMat;
% zerosMat=zeros(size(result));
% zerosMat(resMinusMax==0)=1;
% 
% classOut=nan(size(data.LDR));
% 
% for ii=1:6
%     testMat=squeeze(zerosMat(ii,:,:));
%     classOut(isnan(classOut) & testMat==1)=ii;
% end

if plotIn.plotMax
    plotResMax(data,result,maxAll,plotIn);
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
        
    % No small frozen in warm region
    replaceMat=zeros(size(classOut));
    replaceMat(data.MELTING_LAYER==10 & classOut==9)=1;
        
    % Replace with closest warm pixel
    [oldR oldC]=find(~isnan(classOut) & replaceMat==0 & data.MELTING_LAYER==10);
    [addR addC]=find(replaceMat==1);
    idx = knnsearch([oldR oldC], [addR addC]);
    nearest_OldValue = classOut(sub2ind(size(classOut), oldR(idx), oldC(idx)));
    classOut(sub2ind(size(classOut), addR, addC))=nearest_OldValue;
    
end

end