function [layerAltsOut,layerIndsOut]=zeroDegIso(data)
% Find zero degree isotherms
oddAngles=find(data.elevation>-70 & data.elevation<70);

data.TEMP(:,oddAngles)=nan;
signTemp=sign(data.TEMP);
zeroDeg=diff(signTemp,1);

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(zeroDeg~=0)=1;

zeroSum=sum(zeroDeg,1);
tempDataY=find(zeroSum~=0);
zeroDeg=zeroDeg(:,tempDataY);

zeroAlt=nan(size(zeroDeg));
data.asl=data.asl(:,tempDataY);
zeroAlt(zeroDeg==1)=data.asl(zeroDeg==1);

%% Connect zero degree layers

% Find how many melting layers there are and connect the right ones
zeroTemp=zeroDeg;

numZero=sum(zeroTemp,1,'omitnan');

% Connect layers
layerInds=nan(1,length(tempDataY));
layerAlts=nan(1,length(tempDataY));

for ii=1:length(tempDataY)
    if numZero==0
        continue
    end
    colInds=find(zeroDeg(:,ii)==1);
    colAltsAll=zeroAlt(:,ii);
    colAlts=colAltsAll(colInds);
    
    if ii>1
        prevAlts=layerAlts(:,ii-1);
    end
    
    % Go through all layers
    for jj=1:numZero(ii)
        % First data points or after nans
        if ii==1 | numZero(ii-1)==0
            if ii==1 & jj==1
                layerInds(jj,ii)=colInds(jj);
                layerAlts(jj,ii)=colAlts(jj);
            else
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        else % Find closest altitude
            zeroAltDiffs=abs(prevAlts-colAlts(jj));
            zeroAltDiffs(isnan(zeroAltDiffs))=100000;
            minDiff=min(zeroAltDiffs);
            if minDiff<50 | (numZero(ii)==1 & numZero(ii-1)==1)
                diffInd=find(zeroAltDiffs==minDiff);
                layerInds(diffInd,ii)=colInds(jj);
                layerAlts(diffInd,ii)=colAlts(jj);
            elseif minDiff~=100000
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        end
    end
end

layerAltsOut=nan(size(layerAlts,1),length(data.time));
layerAltsOut(:,tempDataY)=layerAlts;

layerIndsOut=nan(size(layerInds,1),length(data.time));
layerIndsOut(:,tempDataY)=layerInds;

% %% Add adjustment
% zeroAdjustGates=round(zeroAdjustMeters/oneGate);
% 
% layerIndsAdj=nan(size(layerInds));
% layerAltsAdj=nan(size(layerAlts));
% 
% data.elevation=data.elevation(tempDataY);
% 
% layerIndsAdj(:,data.elevation<0)=layerInds(:,data.elevation<0)-zeroAdjustGates;
% layerIndsAdj(:,data.elevation>=0)=layerInds(:,data.elevation>=0)+zeroAdjustGates;
% 
% layerAltsAdj=layerAlts+zeroAdjustMeters;

end