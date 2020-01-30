clear VEL_corr*
clear VEL_ground*
VEL(VEL<-900)=nan;
VEL_RAW(VEL_RAW<-900)=nan;

rad=pi/180;
noise_floor=-107.75;

ground_speed=sqrt(eastward_velocity.^2 + northward_velocity.^2);

[m n]=size(DBZ);

%VEL=VEL_RAW;
VEL_corr=VEL_RAW;
%VEL_corr(DBMVC<noise_floor)=nan;

% Correction added to pointing angle

%Saangria
% tau_suba=3.45;% tilt angle wrt aircraft
% theta_suba=azimuth-1.5;

%IDEAS IV
% tau_suba=tilt -0.12;% tilt angle wrt aircraft
% theta_suba=rotation - .02; % rotation angle wrt aircraft;
% tau_suba=tilt -0.2;% tilt angle wrt aircraft
% theta_suba=rotation - 0.6; % rotation angle wrt aircraft;

% HCR-Test 2014
tau_suba=tilt - 0.;% tilt angle wrt aircraft
theta_suba=rotation - .0; % rotation angle wrt aircraft;

% tau_suba=tilt - 0.0;% tilt angle wrt aircraft
% theta_suba=rotation - .0; % rotation angle wrt aircraft;

%vel_pitch=ground_speed.*sin(pitch2*pi/180);

% compute track direction
clear track
for i=1:m
%     if eastward_velocity(i)>0 && northward_velocity(i)<0
%         track(i)=atan(-northward_velocity(i)/eastward_velocity(i));
%         track(i)=track(i)*180/pi+90;
%     elseif eastward_velocity(i)>0 && northward_velocity(i)>0
%         track(i)=atan(eastward_velocity(i)/northward_velocity(i));
%         track(i)=track(i)*180/pi;
%     elseif eastward_velocity(i)<0 && northward_velocity(i)>0
%         track(i)=atan(northward_velocity(i)/-eastward_velocity(i));
%         track(i)=track(i)*180/pi+270;
%     elseif eastward_velocity(i)<0 && northward_velocity(i)<0
%         track(i)=atan(-northward_velocity(i)/-eastward_velocity(i));
%         track(i)=track(i)*180/pi+180;
%     end
    track(i)=-atan2(northward_velocity(i),eastward_velocity(i))*180/pi +90;
end

heading(heading<0)=heading(heading<0)+360;
track(track<0)=track(track<0)+360;

drift2=track-heading;
drift2(drift<-300)=drift(drift<-300)+360;

% convert to radians
drift_rad=drift*rad;
pitch_rad=pitch*rad;
roll_rad=roll*rad;
theta_suba_rad=theta_suba*rad;
tau_suba_rad=tau_suba*rad;

% Compute y_t following equation 9 Lee et al (1994)
y_subt=-cos(theta_suba_rad+roll_rad).*cos(drift_rad).*cos(tau_suba_rad).*sin(pitch_rad)...
    +sin(drift_rad).*sin(theta_suba_rad+roll_rad).*cos(tau_suba_rad)...
    +cos(pitch_rad).*cos(drift_rad).*sin(tau_suba_rad);

% Compute z following equation 9 Lee et al (1994)
z=cos(pitch_rad).*cos(tau_suba_rad).*cos(theta_suba_rad+roll_rad)+sin(pitch_rad).*sin(tau_suba_rad);

% compute tau_t following equation 11 Lee et al (1994)
tau_subt=asin(y_subt);

% Compute phi following equation 17 Lee et al (1994)
phi=asin(z);

% Compute platform motion based on Eq 27 from Lee et al (1994)
vr_platform=-ground_speed.*sin(tau_subt)-vertical_velocity.*sin(phi);
% vr_platform=-ground_speed.*sin(tau_subt);%-vertical_velocity.*sin(phi);
test=vr_platform;

for i=1:m
for j=1:n
    VEL_corr(i,j)=VEL_corr(i,j)-vr_platform(i);
    while VEL_corr(i,j)<-nyquist_velocity(i)
        VEL_corr(i,j)=VEL_corr(i,j)+2*nyquist_velocity(i);
    end
    while VEL_corr(i,j)>nyquist_velocity(i)
        VEL_corr(i,j)=VEL_corr(i,j)-2*nyquist_velocity(i);
    end

end
end

% Compute the velocity of the ground
DBZ_tmp=DBZ;
DBZ_tmp(:,1:15)=0;
[bla ground_index]=max(DBZ_tmp,[],2);clear DBZ_tmp
% ground_index=ground_index+1; % test

SW=WIDTH;
SW(SW<-1)=nan;
SW(SNRVC<0)=nan;
SW_mean=nanmean(SW,2);

SW_diff2=diff(SW,[],1);
SW_diff2(m,1:n)=0;
SW_diff=nanmean(SW_diff2,2);

VR=VEL;
VR(VR<-15)=nan;
VR(SNRVC<0)=nan;
VR_diff2=diff(VR,[],1);
VR_diff2(m,1:n)=0;
VR_diff=nanmean(VR_diff2,2);

VEL_corr_ground=VEL_corr;
VEL_corr_ground2=VEL;
for i=1:m
    VEL_ground_ucorr(i)=VEL(i,ground_index(i));
    VEL_ground(i)=VEL_corr(i,ground_index(i));
    if abs(VR_diff(i))>0.0;
        VEL_corr_ground(i,:)=VEL_corr(i,:)-VEL_ground(i);
        VEL_corr_ground2(i,:)=VEL(i,:)-VEL_ground_ucorr(i);
    end
end

figure;plot(time,VEL_ground)
hold on;plot(time,VEL_ground_ucorr,'r') %-vertical_velocity,'r')
hold on;plot(time,vr_platform,'k')
%axis([0 70 -1 1])

i_level=find(abs(roll<5) & elevation<-80 & VEL_ground<min(2*std(VEL_ground),1.25));
vr_bias=mean(VEL_ground(i_level))


% num_filter=4;
% [m n]=size(VEL);

% VEL_ground_smooth=VEL_ground;
% for i=1+num_filter:m-num_filter
%     VEL_ground_smooth(i)=mean(VEL_ground(i-num_filter:i+num_filter));
% end

%hold on;plot(time,VEL_ground_smooth,'r')

% for i=1:m
% for j=1:n
%     vel_display2(i,j)=vel_display(i,j)-VEL_ground_smooth(i);
% end
% end
