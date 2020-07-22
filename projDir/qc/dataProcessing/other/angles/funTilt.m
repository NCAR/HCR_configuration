function outFun = funTilt(tilt,roll,pitch,elev,heading,az)
% Function for tilt solver
z=sin(elev);

crotation=(z-sin(pitch).*sin(tilt))./(cos(pitch).*cos(tilt));
srotation=real(sqrt(1-crotation.^2));

x=-crotation.*sin(heading).*cos(tilt).*sin(pitch) + ...
    cos(heading).*srotation.*cos(tilt) + ...
    sin(heading).*cos(pitch).*sin(tilt);

y=-crotation.*cos(heading).*cos(tilt).*sin(pitch) - ...
    sin(heading).*srotation.*cos(tilt) + ...
    cos(pitch).*cos(heading).*sin(tilt);

outFun=x-y.*tan(az);
end

