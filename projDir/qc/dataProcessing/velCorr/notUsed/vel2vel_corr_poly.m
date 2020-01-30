function [VEL_CORR,velSmooth] = vel2vel_corr_poly(data,ptpIn,polyOrder)
% This function calculates VEL_CORR from VEL (it doesn't do the correction
% for the motion, i.e. VEL from VEL_RAW)
%
% Input variables:
% VEL: motion corrected velocity
% DBZ
% range
% alt
%
% Output variables
% VEL_CORR
% VEL_CORR_SURF: corrected velocity at the ground
% rangeToSurf: range from instrument to surface (should be very close to
% alt)

% Check if input dimensions are correct, otherwise flip
flipdim=0;
if size(data.dbz,2)~=length(data.alt)
    data.vel=data.vel';
    data.dbz=data.dbz';
    data.range=data.range';
    flipdim=1;
end

ptp=ptpIn*10/2;

data.vel(data.vel<-900)=nan;

% Find surface indices
data.dbz(1:15,:)=0;

rangeTemp=data.range;
DBZmask=data.dbz; %Mask with DBZ values that are close to the surface

rangeTemp(find(abs(data.range-data.alt)>150))=nan;
DBZmask(isnan(rangeTemp))=nan;
[bla ground_index]=nanmax(DBZmask,[],1);
wrong_ground_ind=find(ground_index==1);

% Convert to linear indices
linInd=sub2ind(size(data.dbz),ground_index,1:length(data.alt));

VEL_ground=data.vel(linInd); %Velocity at the ground
VEL_ground(wrong_ground_ind)=nan;

velSmooth=nan(length(linInd),1);

for ii=(ptp+1):length(data.time)-ptp
    dataInd=ii-ptp:1:ii+ptp;
    grabInd=ptp+1;
    
    timeNum=datenum(data.time(dataInd));
    
    [polFit S Mu]= polyfit(timeNum,VEL_ground(dataInd),polyOrder);
    
    velSmooth(ii)=polyval(polFit,timeNum(grabInd),[],Mu);
    %polVal = polyval(polFit,timeNum,[],Mu);
    
%     close all
%     
%     f1=figure('DefaultAxesFontSize',14);
%     set(f1,'Position',[200 500 1500 500]);
%     hold on
%     plot(data.time,VEL_ground);
%     plot(data.time(dataInd),VEL_ground(dataInd));
%     plot(data.time(dataInd),polVal);
end


VEL_CORR=data.vel-velSmooth';

% Flip output dimensions if necessary
if flipdim
    VEL_CORR=VEL_CORR';
end
end

