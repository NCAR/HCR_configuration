function [velSmoothP,surfInd]=velCorrNadir(data,polyTimePeriod,polyOrder)
%% Fill in extinct echo

disp(['Filling extinct areas.']);
% Remove up pointing
dbzDown=data.DBZ;
dbzDown(:,data.elevation>0)=nan;

velDown=data.VEL;
velDown(:,data.elevation>0)=nan;

% Get surface indices
[surfInd rowInd rangeToSurf] = hcrSurfInds(data);

surfDBZ=dbzDown(surfInd);
surfDBZlin=10.^(surfDBZ./10);

% SurfVel is the surface velocity that is used for the polinomial
% fit
surfVelOrig=velDown(surfInd)';
surfVel=velDown(surfInd);

% We remove all the data that we don't want to include in the fit

% Remove all data where we don't see the surface
surfVel(isnan(surfDBZlin))=nan;

% Remove data that is out of range
surfVel(isnan(rowInd))=nan;

% Calculate standard deviation
surfStd=movstd(surfVel,100,'includenan'); % Sets to nan when there is at least one nan

% Enlarge missing data areas
surfTemp=movmean(surfVel,100,'includenan'); % Sets to nan when there is at least one nan
surfVel(isnan(surfTemp))=nan;

% standard deviation is too high
surfVel(surfDBZlin<10000 & surfStd>0.5)=nan;

% Fill in the missing data with moving average from before the gap
surfMean=movmedian(surfVel,100,'omitnan'); % Ignores nans and uses the rest of the data
surfMean(isnan(surfVel))=nan;

aboveGround=data.altitude-data.TOPO;

surfNan=find(isnan(surfVel) & data.elevation<=0 & aboveGround>10);

for ll=1:length(surfNan)
    if ~isnan(surfMean(surfNan(ll)-1)) % At the beginning of the gap we have moving average data
        surfVel(surfNan(ll))=surfMean(surfNan(ll)-1);
    else % Once the moving average turns nan we just keep the previous one going
        surfVel(surfNan(ll))=surfVel(surfNan(ll)-1);
    end
end
%% Make poly fit

disp('Making fit.');

velSmoothPorig=vel2vel_corr_testPoly(surfVel,data.time,polyTimePeriod,polyOrder);

% Fill in data where we don't have polyfit data
velSmoothP=movmedian(surfVelOrig,100,'omitnan');
velSmoothP(~isnan(velSmoothPorig))=velSmoothPorig(~isnan(velSmoothPorig));

end