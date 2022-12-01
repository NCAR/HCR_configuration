function cloudClass=findCloudClass(jEcho,cloudID,jTemp,jMelt,jElev,jTopo,jAsl)
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

% Set elev to strat
% Minimum altitude for low/mid boundary is 2km
% Minimum altitude for mid/high boundary is 4km
distAslTopo=jAsl-jTopo;
melt=jMelt;
melt(distAslTopo<2000)=10;
melt(isnan(jMelt))=nan;
jTemp(distAslTopo<4000 & jTemp<-25)=-25;
% Low
jEcho(jEcho==32 & melt<20)=14;
% Mid
jEcho(jEcho==32 & melt>=20 & jTemp>=-25)=16;
% High
jEcho(jEcho==32 & melt>=20 & jTemp<-25)=18;

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
        % Check for low altitude flight legs
%         altMat=distAslTopo(min(r):max(r),min(c):max(c));
%         altMat(idMat~=ii)=nan;

        cloudMatAlt=cloudMat(1,:);
        cloudMatAlt(distAslTopo(18,min(c):max(c))<500 & jElev(min(c):max(c))>0)=nan;

        inCloudFrac=sum(~isnan(cloudMatAlt))./size(cloudMat,2);
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

    highPix=length(find(cloudMat==18));
    midPix=length(find(cloudMat==16));
    lowPix=length(find(cloudMat==14));

    % Check if convective
    % Calculate convective fraction
    convPix=sum(sum(cloudMat>32));
    stratPix=sum(sum(cloudMat<20));

    convFrac=convPix/(stratPix+convPix);

    if convFrac>0.05 | max(cloudMat(:))==38
        % Young
        if convFrac>0.5
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
                if highPix/length(cloudMat)>0.3
                    cloudClass(linInds)=0;
                else
                    cloudClass(linInds)=32; % Mid
                end
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==34
                if highPix>0 | midPix/length(cloudMat)>0.5
                    cloudClass(linInds)=0;
                else
                    cloudClass(linInds)=31; % Shallow
                end
            end
        end
        %% Stratiform
    else
        % Set the few convective pixels to low numbers
        cloudMat(cloudMat>=20)=0;
        %% Check if precipitating
        topoDistMat=distAslTopo(min(r):max(r),min(c):max(c));
        topoDistMat(isnan(cloudMat))=nan;
        minDist=min(reshape(topoDistMat,1,[]),[],'omitnan');
        
        % Precipitating
        if minDist<500
            if max(reshape(cloudMat,1,[]),[],'omitnan')==18 
                cloudClass(linInds)=13; % Deep
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==16
                cloudClass(linInds)=12; % Mid
            elseif max(reshape(cloudMat,1,[]),[],'omitnan')==14
                cloudClass(linInds)=11; % Shallow
            end

        % Non precipitating
        else
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