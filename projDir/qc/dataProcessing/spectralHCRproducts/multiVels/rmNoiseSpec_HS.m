function powerRMnoiseOut=rmNoiseSpec_HS(powerAdj,vNoise)
% Find mean noise and noise threshold with following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2
powerRMnoiseOut=nan(size(powerAdj));

sampleNumR=size(powerAdj,2);

% Moving average
meanOverPoints=55; % Average over this number of points
movAv=movmedian(powerAdj,meanOverPoints,2);

% pointInds=1:meanOverPoints:sampleNumR;
% pointAv=movAv(:,pointInds);

startNoise=10^((vNoise+10)/10);

powerRMnoiseAv=nan(size(powerAdj));

loopInds=find(any(~isnan(powerAdj),2));
for aa=1:size(loopInds,1)
    ii=loopInds(aa);

    thisMov=movAv(ii,:);

    % Find noise floor and noise threshold
    thisMovLin=10.^(thisMov./10);
    [noiseThreshM,meanNoiseM]=findNoiseThresh(thisMovLin,5,startNoise);

    % Remove noise
    thisMovRM=thisMov;
    thisMovRM(thisMovRM<noiseThreshM)=nan;
    thisMask=~isnan(thisMovRM);
    thisMask=bwareafilt(thisMask,1);
    thisMovRM(thisMask==0)=nan;
    powerRMnoiseAv(ii,:)=thisMovRM;
           
    plot(1:sampleNumR,powerAdj(ii,:),'-b')
    hold on
    plot(1:sampleNumR,thisMov,'-g','LineWidth',1.5);
    plot(1:sampleNumR,thisMovRM,'-r','LineWidth',1.5);
    plot([1,sampleNumR],[noiseThreshM,noiseThreshM],'-c','LineWidth',2);
    plot([1,sampleNumR],[meanNoiseM,meanNoiseM],'-k','LineWidth',2);
    xlim([1,sampleNumR]);
    hold off

    if ii==97
        stop1=1;
    end

end
end