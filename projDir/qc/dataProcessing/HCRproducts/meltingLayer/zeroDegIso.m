function [layerAltsOut,layerIndsOut,layerVelsOut,meltOnly,tempOut]=zeroDegIso(data)
% Find zero degree isotherms
oddAngles=find(data.elevation>-70 & data.elevation<70);

data.TEMP(:,oddAngles)=nan;
tempOut=data.TEMP>0;
%tempOut=imclose(tempOut,strel('disk',3));
tempOut=double(tempOut);
tempOut(tempOut==0)=-1;
tempOut(isnan(data.TEMP))=nan;

signTemp=sign(tempOut);
zeroDeg=diff(signTemp,1);

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(:,data.elevation>0)=-zeroDeg(:,data.elevation>0);
zeroDegUpDown=zeroDeg;
zeroDeg(zeroDeg~=0)=1;

zeroSum=sum(zeroDeg,1);
tempDataY=find(zeroSum~=0);
zeroDeg=zeroDeg(:,tempDataY);
zeroDegUpDown=zeroDegUpDown(:,tempDataY);

zeroAlt=nan(size(zeroDeg));
data.asl=data.asl(:,tempDataY);
zeroAlt(zeroDeg==1)=data.asl(zeroDeg==1);

zeroVel=nan(size(zeroDeg));
data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);
data.VEL_MASKED=data.VEL_MASKED(:,tempDataY);
zeroVel(zeroDeg==1)=data.VEL_MASKED(zeroDeg==1);

%% Connect zero degree layers

% Find how many melting layers there are and connect the right ones
zeroTemp=zeroDeg;

numZero=sum(zeroTemp,1,'omitnan');

% Connect layers
layerInds=nan(1,length(tempDataY));
layerAlts=nan(1,length(tempDataY));
layerVels=nan(1,length(tempDataY));
meltOnly=[];

for ii=1:length(tempDataY)
    if numZero==0
        continue
    end
    colInds=find(zeroDeg(:,ii)==1);
    zUD=zeroDegUpDown(colInds,ii);
    colAltsAll=zeroAlt(:,ii);
    colAlts=colAltsAll(colInds);
    colVelsAll=zeroVel(:,ii);
    colVels=colVelsAll(colInds);
    
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
                layerVels(jj,ii)=colVels(jj);
                if  zUD(jj)>0
                    meltOnly=[meltOnly;1];
                else
                    meltOnly=[meltOnly;0];
                end
            else
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerVels=cat(1,layerVels,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
                layerVels(size(layerInds,1),ii)=colVels(jj);
                if  zUD(jj)>0
                    meltOnly=[meltOnly;1];
                else
                    meltOnly=[meltOnly;0];
                end
            end
        else % Find closest altitude
            zeroAltDiffs=abs(prevAlts-colAlts(jj));
            zeroAltDiffs(isnan(zeroAltDiffs))=100000;
            minDiff=min(zeroAltDiffs);
            if minDiff<50 | (numZero(ii)==1 & numZero(ii-1)==1)
                diffInd=find(zeroAltDiffs==minDiff);
                layerInds(diffInd,ii)=colInds(jj);
                layerAlts(diffInd,ii)=colAlts(jj);
                layerVels(diffInd,ii)=colVels(jj);
            elseif minDiff~=100000
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerVels=cat(1,layerVels,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
                layerVels(size(layerInds,1),ii)=colVels(jj);
                if  zUD(jj)>0
                    meltOnly=[meltOnly;1];
                else
                    meltOnly=[meltOnly;0];
                end
            end
        end
    end
end

layerAltsOut=nan(size(layerAlts,1),length(data.time));
layerAltsOut(:,tempDataY)=layerAlts;

layerVelsOut=nan(size(layerVels,1),length(data.time));
layerVelsOut(:,tempDataY)=layerVels;

layerIndsOut=nan(size(layerInds,1),length(data.time));
layerIndsOut(:,tempDataY)=layerInds;

end