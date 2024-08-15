function momentsSpec=calcMomentsSpec_higherMoments(specPowerDB,specPhase,ii,momentsSpec)

x=specPhase;
y=10.^(specPowerDB./10);

% VEL
momentsSpec.velRaw(:,ii)=sum(y.*x,2,'omitnan')./sum(y,2,'omitnan');

% WIDTH
momentsSpec.width(:,ii)=real((sum(y.*(x-momentsSpec.velRaw(:,ii)).^2,2,'omitnan')./sum(y,2,'omitnan')).^0.5);

% SKEWNESS
momentsSpec.skew(:,ii)=sum(y.*(x-momentsSpec.velRaw(:,ii)).^3,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^3);

% KURTOSIS
momentsSpec.kurt(:,ii)=(sum(y.*(x-momentsSpec.velRaw(:,ii)).^4,2,'omitnan')./(sum(y,2,'omitnan').*momentsSpec.width(:,ii).^4))-3;

end