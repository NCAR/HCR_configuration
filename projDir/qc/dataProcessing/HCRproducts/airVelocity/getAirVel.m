function airVel=getAirVel(powerAdj,phaseAdj,sampleNum)
% Get air velocity from spectral data

airVel=nan(size(powerAdj,1),1);

powerMed=movmedian(powerAdj,round(sampleNum/10),2);

for ii=1:length(airVel)

    plot(phaseAdj(ii,:),powerAdj(ii,:));
    hold on
    plot(phaseAdj(ii,:),powerMed(ii,:));
    hold off
    ylim([-100,-40]);

end