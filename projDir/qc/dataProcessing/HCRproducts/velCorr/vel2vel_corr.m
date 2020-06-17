function [VEL_CORR,VEL_CORR_SURF,rangeToSurf,velSmoothOut] = vel2vel_corr(VEL,DBZ,range,alt)
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
if size(DBZ,2)~=length(alt)
    VEL=VEL';
    DBZ=DBZ';
    range=range';
    flipdim=1;
end

VEL(VEL<-900)=nan;

% Find surface indices
DBZ(1:15,:)=0;

rangeTemp=range;
DBZmask=DBZ; %Mask with DBZ values that are close to the surface

rangeTemp(find(abs(range-alt)>150))=nan;
DBZmask(isnan(rangeTemp))=nan;
[bla ground_index]=nanmax(DBZmask,[],1);
wrong_ground_ind=find(ground_index==1);

% Convert to linear indices
linInd=sub2ind(size(DBZ),ground_index,1:length(alt));

rangeToSurf=range(linInd); %Range to surface, should be very close to alt
rangeToSurf(wrong_ground_ind)=nan;
VEL_ground=VEL(linInd); %Velocity at the ground
VEL_ground(wrong_ground_ind)=nan;

% First fir filter, remove spikes
num_iterations1=10;
VG_fir=nan(num_iterations1,length(ground_index));
VG_fir(1,:)=Hubbert_filter(VEL_ground);
for i=2:num_iterations1
    VG_fir(i,:)=Hubbert_filter(VG_fir(i-1,:));
end

VG_fir_tofilter=VG_fir(num_iterations1,:);
VG_fir_inds=find(abs(VEL_ground-VG_fir(num_iterations1,:))<0.2);
VG_fir_tofilter(VG_fir_inds)=VEL_ground(VG_fir_inds);

% Second fir filter
num_iterations2=3;
velSmooth=nan(num_iterations2,length(ground_index));
velSmooth(1,:)=Hubbert_filter(VG_fir_tofilter);
for i=2:num_iterations2
    velSmooth(i,:)=Hubbert_filter(velSmooth(i-1,:));
end

VEL_CORR=VEL-velSmooth(num_iterations2,:);
VEL_CORR_SURF=VEL_CORR(linInd);
VEL_CORR_SURF(wrong_ground_ind)=nan;

velSmoothOut=velSmooth(num_iterations2,:);

% Flip output dimensions if necessary
if flipdim
    VEL_CORR=VEL_CORR';
end
end

