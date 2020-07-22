function [azEl_cfrad,azEl_lee,azEl_field] = azEl(roll,pitch,heading,rotation,tilt,drift)
%Calculate azimuth and elevation angle.
xyz_cfrad=[];
xyz_lee=[];
xyz_field=[];

roll=deg2rad(roll);
pitch=deg2rad(pitch);
heading=deg2rad(heading);
rotation=deg2rad(rotation);
tilt=deg2rad(tilt);
drift=deg2rad(drift);

for ii=1:length(roll)
    % From the equation in the cfradial doc
    Mr=[cos(roll(ii)),0,sin(roll(ii));...
        0,1,0;...
        -sin(roll(ii)),0,cos(roll(ii))];
    
    Mp=[1,0,0;...
        0,cos(pitch(ii)),-sin(pitch(ii));...
        0,sin(pitch(ii)),cos(pitch(ii))];
    
    Mh=[cos(heading(ii)),sin(heading(ii)),0;...
        -sin(heading(ii)),cos(heading(ii)),0;...
        0,0,1];
    
    MhMpMr=Mh*Mp*Mr;
    
    angvec_plane=[sin(rotation(ii))*cos(tilt(ii));...
        sin(tilt(ii));...
        cos(rotation(ii))*cos(tilt(ii))];
    
    xyzdivr=MhMpMr*angvec_plane;
    
    xyz_cfrad=cat(1,xyz_cfrad,xyzdivr');
    
    % To compare with Lee et al:
    xyzdivr_2=[-cos(rotation(ii)+roll(ii))*sin(heading(ii))*cos(tilt(ii))*sin(pitch(ii)) + ...
        cos(heading(ii))*sin(rotation(ii)+roll(ii))*cos(tilt(ii)) + ...
        sin(heading(ii))*cos(pitch(ii))*sin(tilt(ii)); ...
        -cos(rotation(ii)+roll(ii))*cos(heading(ii))*cos(tilt(ii))*sin(pitch(ii)) - ...
        sin(heading(ii))*sin(rotation(ii)+roll(ii))*cos(tilt(ii)) + ...
        cos(pitch(ii))*cos(heading(ii))*sin(tilt(ii));...
        cos(pitch(ii))*cos(tilt(ii))*cos(rotation(ii)+roll(ii)) + ...
        sin(pitch(ii))*sin(tilt(ii))];
    
    xyz_lee=cat(1,xyz_lee,xyzdivr_2');
    
    % To compare with the code that is actually used
    
    xsubt=cos(rotation(ii)+roll(ii)) * sin(drift(ii)) * cos(tilt(ii)) * sin(pitch(ii)) + ...
        cos(drift(ii)) * sin(rotation(ii)+roll(ii)) * cos(tilt(ii)) - ...
        sin(drift(ii)) * cos(pitch(ii)) * sin(tilt(ii));
    
    ysubt = -cos(rotation(ii)+roll(ii)) * cos(drift(ii)) * cos(tilt(ii)) * sin(pitch(ii)) + ...
        sin(drift(ii)) * sin(rotation(ii)+roll(ii)) * cos(tilt(ii)) + ...
        cos(pitch(ii)) * cos(drift(ii)) * sin(tilt(ii));
    
    zsubt = cos(pitch(ii)) * cos(tilt(ii))* cos(rotation(ii)+roll(ii)) + ...
        sin(pitch(ii)) * sin(tilt(ii));
    
    xyzdivr_3=[xsubt;ysubt;zsubt];
    
    xyz_field=cat(1,xyz_field,xyzdivr_3');
    
end

%cfradial doc

az_cfrad_rad=atan2(xyz_cfrad(:,1),xyz_cfrad(:,2));
az_cfrad=rad2deg(mod(az_cfrad_rad,2*pi));

el_cfrad=rad2deg(asin(xyz_cfrad(:,3)));

azEl_cfrad=cat(2,az_cfrad,el_cfrad);

% Lee et al.

az_lee_rad=atan2(xyz_lee(:,1),xyz_lee(:,2));
az_lee=rad2deg(mod(az_lee_rad,2*pi));

el_lee=rad2deg(asin(xyz_lee(:,3)));

azEl_lee=cat(2,az_lee,el_lee);

% Field

az_t=atan2(xyz_field(:,1),xyz_field(:,2));
el_field=rad2deg(asin(xyz_field(:,3)));
    
track=heading+drift;
az_field=rad2deg(mod(az_t+track,2*pi));

azEl_field=cat(2,az_field,el_field);

end

