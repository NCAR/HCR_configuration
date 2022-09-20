function cloudClass=findCloudClass(jEcho,cloudID,jTemp,jElev,jTopo,jAsl)
% Classify clouds based on their conv/strat properties

% Not classified=0
% stratLow=1
% stratMid=2;
% stratHigh=3
% stratPrecipShallow=11
% stratPrecipMid=12
% stratPrecipDeep=13
% convYoungShallow=21
% convYoungMid=22
% convYongDeep=23
% convMatureShallow=31
% convMatureMid=32
% convMatureDeep=33

cloudClass=nan(size(jEcho));

% Set elev to mixed
jEcho(jEcho==32)=25;

% Calculate alt minus topo
topoDist=jAsl-jTopo;

% Loop through clouds
for ii=1:max(reshape(cloudID,1,[]))

    linInds=find(cloudID==ii);
    [r,c]=ind2sub(size(cloudID),linInds);

    idMat=cloudID(min(r):max(r),min(c):max(c));
    cloudMat=jEcho(min(r):max(r),min(c):max(c));
    cloudMat(idMat~=ii)=nan;

    %% Sort out clouds we don't classify

    % If only convective with no sub class, don't use
    if any(cloudMat==30)
        cloudClass(linInds)=0;
        continue
    end

    % If mostly mixed, don't use
    mixedPix=length(find(cloudMat==25));
    mixedFrac=mixedPix./(sum(sum(~isnan(cloudMat))));

    if mixedFrac>0.8
        cloudClass(linInds)=0;
        continue
    end

    cloudMat(cloudMat==25)=0;
    %% Calculate in-cloud fraction

    inCloudFrac=0;

    if (min(r))==18
        inCloudFrac=sum(~isnan(cloudMat(1,:)))./size(cloudMat,2);
        if inCloudFrac>0.3
            % Pointing up
            if median(jElev(min(c):max(c)))>0
                cloudClass(linInds)=0;
                continue
                % Check temperature
            elseif mean(jTemp(min(r),min(c):max(c)),'omitnan')>-30
                cloudClass(linInds)=0;
                continue
            end
        end
    end

    %% Convective
    % Check if convective
    if max(reshape(cloudMat,1,[]),[],'omitnan')>32
        % Calculate convective fraction
        convPix=sum(sum(cloudMat>32));
        stratPix=sum(sum(cloudMat<20));

        convFrac=convPix/(stratPix+convPix);

        % Young
        if convFrac>0.7
            if max(reshape(cloudMat,1,[]),[],'omitnan')==38 
                cloudClass(linInds)=23; % Deep
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==36
                cloudClass(linInds)=22; % Mid
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==34
                cloudClass(linInds)=21; % Shallow
            end
        else % Mature
            if max(reshape(cloudMat,1,[]),[],'omitnan')==38
                cloudClass(linInds)=33; % Deep
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==36
                cloudClass(linInds)=32; % Mid
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==34
                cloudClass(linInds)=31; % Shallow
            end
        end
        %% Stratiform
    else
        %% Check if precipitating
        topoDistMat=topoDist(min(r):max(r),min(c):max(c));
        topoDistMat(isnan(cloudMat))=nan;
        minDist=min(reshape(topoDistMat,1,[]),[],'omitnan');
        
        % Precipitating
        if minDist<200
            if max(reshape(cloudMat,1,[]),[],'omitnan')==18 
                cloudClass(linInds)=13; % Deep
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==16
                cloudClass(linInds)=12; % Mid
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==14
                cloudClass(linInds)=11; % Shallow
            end

        % Non precipitating
        else
            highPix=length(find(cloudMat==18));
            midPix=length(find(cloudMat==16));
            lowPix=length(find(cloudMat==14));

            [maxPix,maxInd]=max([highPix,midPix,lowPix]);

            if maxInd==1
                cloudClass(linInds)=3; % High
            elseif maxInd==2
                cloudClass(linInds)=2; % Mid
            else
                cloudClass(linInds)=1; % Low
            end
        end
    end

end