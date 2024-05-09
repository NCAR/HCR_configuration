function momentsSpec=calcMomentsSpec_higherMoments(specPowerDB,specPhase,ii,momentsSpec,data)

noiseLinV=10.^(data.noise_v./10);
x=specPhase;
y=10.^(specPowerDB./10)-noiseLinV;

% VEL
momentsSpec.velRaw(:,ii)=sum(y.*x,2,'omitnan')./sum(y,2,'omitnan');

% Correct velocity for aircraft motion
% xCorr=sind(momentsSpec.azimuth_vc(ii)).*cosd(momentsSpec.elevation(ii)).*momentsSpec.eastward_velocity(ii);
% yCorr=cosd(momentsSpec.azimuth_vc(ii)).*cosd(momentsSpec.elevation(ii)).*momentsSpec.northward_velocity(ii);
% zCorr=sind(momentsSpec.elevation(ii)).*momentsSpec.vertical_velocity(ii);
% momentsSpec.vel(:,ii)=momentsSpec.velRaw(:,ii)+xCorr+yCorr+zCorr;

% WIDTH
momentsSpec.width(:,ii)=real((sum(y.*(x-momentsSpec.velRaw(:,ii)).^2,2,'omitnan')./sum(y,2,'omitnan')).^0.5);

% % Correct width for aircraft motion
% velAircraft=sqrt(momentsSpec.eastward_velocity(ii).^2+momentsSpec.northward_velocity(ii).^2);
% deltaC=0.3.*velAircraft.*sin(deg2rad(momentsSpec.elevation(ii))).*deg2rad(data.beamwidth_v);
% momentsSpec.width(:,ii)=real(sqrt(widthRaw.^2-deltaC.^2));

% SKEWNESS
momentsSpec.skew(:,ii)=sum(y.*(x-momentsSpec.velRaw(:,ii)).^3,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^3);

% KURTOSIS
momentsSpec.kurt(:,ii)=(sum(y.*(x-momentsSpec.velRaw(:,ii)).^4,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^4))-3;

end