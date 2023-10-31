function momentsSpec=calcMomentsSpec_higherMoments(specPowerDB,specPhase,ii,momentsSpec,data)

noiseLinV=10.^(data.noise_v./10);
x=specPhase;
y=10.^(specPowerDB./10)-noiseLinV;

% VEL
momentsSpec.vel(:,ii)=sum(y.*x,2,'omitnan')./sum(y,2,'omitnan');

% WIDTH
momentsSpec.width(:,ii)=real((sum(y.*(x-momentsSpec.vel(:,ii)).^2,2,'omitnan')./sum(y,2,'omitnan')).^0.5);

% SKEWNESS
momentsSpec.skew(:,ii)=sum(y.*(x-momentsSpec.vel(:,ii)).^3,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^3);

% KURTOSIS
momentsSpec.kurt(:,ii)=sum(y.*(x-momentsSpec.vel(:,ii)).^4,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^4);

end