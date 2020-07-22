function [VEL_ground ground_speed] = calcVelCorr(drift,pitch,roll,rotation,tilt,eastward_velocity,...
    northward_velocity,vertical_velocity,VEL_RAW,DBZ)

% Convert to radian
drift_rad=deg2rad(drift);
pitch_rad=deg2rad(pitch);
roll_rad=deg2rad(roll);
theta_suba_rad=deg2rad(rotation);
tau_suba_rad=deg2rad(tilt);

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
ground_speed=sqrt(eastward_velocity.^2 + northward_velocity.^2);
vr_platform=-ground_speed.*sin(tau_subt)-vertical_velocity.*sin(phi);

VEL_corr=VEL_RAW-vr_platform;

% Compute the velocity of the ground
DBZ_tmp=DBZ;
DBZ_tmp(1:15,:)=0;
[bla ground_index]=max(DBZ_tmp,[],1);

ground_index(ground_index==1)=nan;

VEL_ground_raw=nan(1,size(VEL_RAW,2));
VEL_ground=nan(1,size(VEL_RAW,2));

for jj=1:size(VEL_RAW,2)
    if ~isnan(ground_index(jj))
        if DBZ(ground_index(jj),jj)>20
            VEL_ground_raw(jj)=VEL_RAW(ground_index(jj),jj);
            VEL_ground(jj)=VEL_corr(ground_index(jj),jj);
        end
    end
end
end
