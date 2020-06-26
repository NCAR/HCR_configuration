% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function [BBfinished]= f_meltLayer_altOnly(data)
%% Find zero degree altitude

disp('Searching 0 deg altitude ...');

zeroAdjustMeters=350;

oneGate=data.range(2)-data.range(1);
zeroAdjustGates=round(zeroAdjustMeters/oneGate);

timeMat=repmat(data.time,size(data.TEMP,1),1);

tempTemp=data.TEMP;
tempTemp(1:16,:)=nan;
signTemp=sign(tempTemp);
zeroDeg=diff(signTemp,1);

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(zeroDeg~=0)=1;

zeroAlt=nan(size(data.asl));
zeroAlt(zeroDeg==1)=data.asl(zeroDeg==1);

%% Connect zero degree layers

disp('Connecting zero degree layers ...');

% Find how many melting layers there are and connect the right ones
zeroTemp=zeroDeg;

% % Remove contiguous layers that are too small
% if length(data.time)>9000
%     zeroCut=3000;
%     CC = bwconncomp(zeroTemp);
%
%     for ii=1:CC.NumObjects
%         area=CC.PixelIdxList{ii};
%         if length(area)<=zeroCut
%             zeroTemp(area)=nan;
%         end
%     end
% end

numZero=sum(zeroTemp,1,'omitnan');

% Connect layers
layerInds=nan(1,length(data.time));
layerAlts=nan(1,length(data.time));

for ii=1:length(data.time)
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
                layerInds=cat(1,layerInds,nan(size(data.time)));
                layerAlts=cat(1,layerAlts,nan(size(data.time)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        else % Find closest altitude
            zeroAltDiffs=abs(prevAlts-colAlts(jj));
            zeroAltDiffs(isnan(zeroAltDiffs))=100000;
            minDiff=min(zeroAltDiffs);
            if minDiff<50
                diffInd=find(zeroAltDiffs==minDiff);
                layerInds(diffInd,ii)=colInds(jj);
                layerAlts(diffInd,ii)=colAlts(jj);
            elseif minDiff~=100000;
                % Add new row
                layerInds=cat(1,layerInds,nan(size(data.time)));
                layerAlts=cat(1,layerAlts,nan(size(data.time)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        end
    end
end

%% Add adjusted layers
layerIndsAdj=nan(size(layerInds));
layerAltsAdj=nan(size(layerAlts));

adjustMeters=zeroAdjustGates*oneGate;

for kk=1:length(data.time)
    if data.elevation(kk)<0
        layerIndsAdj(:,kk)=layerInds(:,kk)+zeroAdjustGates;
    else
        layerIndsAdj(:,kk)=layerInds(:,kk)-zeroAdjustGates;
    end
    layerAltsAdj(:,kk)=layerAlts(:,kk)-adjustMeters;
end

%% Remove data that is not suitable
LDRdata=data.LDR;
LDRdata(data.FLAG>1)=nan;
LDRdata(LDRdata<-16 | LDRdata>-7)=nan;
LDRdata(data.range<150)=nan;

VELdata=data.VEL_CORR;
VELdata(data.FLAG>1)=nan;
VELdata(data.range<150)=nan;

%% Tightened backlobe
% Initiate mask
blMask=zeros(size(data.DBZ));

blMask(data.DBZ<-18 & data.WIDTH>1)=1;

% Only within right altitude
rightAlt=data.altitude-data.TOPO;

altMat=repmat(rightAlt,size(data.range,1),1);
% Lower limit
blMask(data.range<(altMat-100))=0;
% Upper limit
blMask(data.range>(altMat+2500))=0;

% Only when scanning up
blMask(:,find(data.elevation<0))=0;

LDRdata(blMask==1)=nan;
VELdata(blMask==1)=nan;

%% Find altitude of bright band

BBaltAll={};
timeIall={};
BBaltInterpAll={};
BBaltZeroAll={};
timeIZeroAll={};

transLength=6000;

for kk=1:size(layerAltsAdj,1)
    timeInds=find(~isnan(layerAltsAdj(kk,:)));
    layerAltsTemp=layerAltsAdj(kk,timeInds);
    rowInds=layerIndsAdj(kk,timeInds);
    
    % LDR
    % Find maximum LDR level in area around zero degree altitude
    maxLevelLDR=nan(size(rowInds));
    for ii=1:length(rowInds)
        vertColLDR=LDRdata(:,timeInds(ii));
        vertColLDR(rowInds(ii)+12:end)=nan;
        vertColLDR(1:max([rowInds(ii)-12,1]))=nan;
        
        % Check if all nan
        if min(isnan(vertColLDR))==0
            %vertColData=vertColLDR(~isnan(vertColLDR));
            maxLevelLDR(ii)=min(find(vertColLDR==nanmax(vertColLDR)));
        end
    end
    
    if min(isnan(maxLevelLDR))==0
        ASLlayer=data.asl(:,timeInds);
           
        % LDR       
        colIndsLDR=1:1:size(ASLlayer,2);
        linearIndLDR = sub2ind(size(ASLlayer), maxLevelLDR, colIndsLDR);
        linearIndLDR(isnan(maxLevelLDR))=[];
        
        % Raw altitude of melting layer
        LDRaltRawIn=ASLlayer(linearIndLDR);
        LDRaltRaw=nan(size(rowInds));
        LDRaltRaw(find(~isnan(maxLevelLDR)))=LDRaltRawIn;
        % Mean altitude of meltig layer
        LDRalt=movmedian(LDRaltRaw,300,'omitnan');
        LDRalt(isnan(LDRaltRaw))=nan;
        % Standard deviation
        LDRaltS=movstd(LDRaltRaw,300,'omitnan');
        LDRaltS(isnan(LDRaltRaw))=nan;
        % Remove data with too much std
        LDRaltRaw(LDRaltS>100)=nan;
        % Distance between raw and mean altitude
        LDRloc=abs(LDRaltRaw-LDRalt);
        
        % Remove data where distance is more than 200 m
        LDRaltRaw(LDRloc>50)=nan;
        
        % Adjust zero degree layer
        zeroDist=LDRalt-layerAltsTemp;
        
        % Remove small data stretches
        zeroMask=zeros(size(zeroDist));
        zeroMask(~isnan(zeroDist))=1;
        
        zeroMask=bwareaopen(zeroMask,5);
        zeroDist(zeroMask==0)=nan;
        
        zeroDistSmooth=movmean(zeroDist,3000,'omitnan');
        zeroDistSmooth(isnan(zeroDist))=nan;        
        
        noNanDist=fillmissing(zeroDistSmooth,'linear','EndValues','nearest');
        if sum(~isnan(noNanDist))>0
            layerAltsTemp=layerAltsTemp+noNanDist;
            layerAltsAdj(kk,timeInds)=layerAltsTemp;
            
            distInds=round(noNanDist./oneGate);
            rowInds(data.elevation(timeInds)>0)=rowInds(data.elevation(timeInds)>0)-distInds(data.elevation(timeInds)>0);
            rowInds(data.elevation(timeInds)<=0)=rowInds(data.elevation(timeInds)<=0)+distInds(data.elevation(timeInds)<=0);
        end
    end
    
    % VEL
    % Find vel melting layer level
    maxLevelVEL=nan(size(rowInds));
    velMasked=VELdata;
    velMasked(:,data.elevation<0)=-velMasked(:,data.elevation<0);
    
    velTime=velMasked(:,timeInds);
    for ii=1:length(rowInds)
        velTime(rowInds(ii)+50:end,ii)=nan;
        velTime(1:max([rowInds(ii)-50,1]),ii)=nan;
    end
    
    % Find discontinuity in vel field
    % Smooth in range direction
    velSmooth=movmedian(velTime,10,1);
    
    % Calculate velocity steps
    [velSteps,S1,S2] = ischange(velSmooth,1,'MaxNumChanges',1);
    
    % Difference between steps
    diffS1=nan(1,size(velSteps,2));
    varS2=nan(2,size(velSteps,2));
    for ii=1:size(velSteps,2)
        uS1=unique(S1(:,ii));
        uS1(isnan(uS1))=[];
        uS2=unique(S2(:,ii));
        uS2(isnan(uS2))=[];
        
        if length(uS1)==2
            diffS1(ii)=uS1(2)-uS1(1);
        end
        if length(uS2)==2
            varS2(:,ii)=uS2;
        end
    end
    
    maxVar=max(varS2,[],1);
    
    % Find altitude of step
    stepInLin=find(velSteps==1);
    [stepInR,stepInC]=ind2sub(size(velSmooth),stepInLin);
    
    maxLevelVEL=nan(1,length(timeInds));
    maxLevelVEL(stepInC)=stepInR; 
    
%     fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);
%     colmap=jet;
%     colormap(flipud(colmap));
%     subplot(3,1,1)
%     hold on
%     surf(velSmooth,'edgecolor','none');
%     view(2)
%     %colorbar
%     caxis([-5 5])
%     plot(maxLevelVEL,'-b')
%     
%     subplot(3,1,2)
%     hold on
%     plot(diffS1)
%     plot(movmean(diffS1,50,'omitnan'),'linewidth',2)
%     
%     subplot(3,1,3)
%     hold on
%     plot(maxVar);
%     plot(movmean(maxVar,50,'omitnan'),'linewidth',2)
    
    % Remove data that doesn't cut it
    if min(isnan(maxLevelLDR))==0 | min(isnan(maxLevelLDR))==0
                
        % VEL
        colIndsVEL=1:1:size(ASLlayer,2);
        linearIndVEL = sub2ind(size(ASLlayer), maxLevelVEL, colIndsVEL);
        linearIndVEL(isnan(maxLevelVEL))=[];
        
        % Raw altitude of melting layer
        VELaltRawIn=ASLlayer(linearIndVEL);
        VELaltRaw=nan(size(rowInds));
        VELaltRaw(find(~isnan(maxLevelVEL)))=VELaltRawIn;
        
        % Compare with zero degree layer
        VELzeroDiff=VELaltRaw-layerAltsTemp;
        VELaltRaw(abs(VELzeroDiff)>200)=nan;
        
        % Mean altitude of meltig layer
        VELalt=movmedian(VELaltRaw,1000,'omitnan');
        VELalt(isnan(VELaltRaw))=nan;
        
        % Distance between raw and mean altitude
        VELloc=abs(VELaltRaw-VELalt);        
        % Remove data where distance is more than 200 m
        VELaltRaw(VELloc>100)=nan;
        
        % Standard deviation
        VELaltS=movstd(VELaltRaw,300,'omitnan');
        VELaltS(isnan(VELaltRaw))=nan;
        % Remove data with too much std
        VELaltRaw(VELaltS>100)=nan;
                
        % Combine max levels
        BBaltRaw=VELaltRaw;
        BBaltRaw(~isnan(LDRaltRaw))=LDRaltRaw(~isnan(LDRaltRaw));        
        
        % Remove data that is too short
        BBmask=zeros(size(BBaltRaw));
        BBmask(~isnan(BBaltRaw))=1;
        
        indMask=find(abs(diff(BBmask))==1);
        diffMask=diff(indMask);
        
        diffOnes=find(diffMask==1);
        indsOnes=indMask(diffOnes);
        
        for ii=1:length(indsOnes)
            if BBmask(indsOnes(ii)+1)==0
                BBmask(indsOnes(ii)+1)=1;
            end
        end
        
        % Keep data that has only a few nan
        zeroCut=20;
        CC = bwconncomp(BBmask);
        
        for ii=1:CC.NumObjects
            area=CC.PixelIdxList{ii};
            if length(area)<=zeroCut
                BBaltRaw(area)=nan;
            end
        end
                
        % Avoid jumps between gates
        BBaltRaw=movmedian(BBaltRaw,20,'omitnan');
        
        % Interpolate between good values
        BBaltInterp=nan(size(BBaltRaw));
        BBaltZero=nan(size(BBaltRaw));
        
        if min(isnan(BBaltRaw))==0
            % Mask
            maskBBalt=zeros(size(BBaltRaw));
            maskBBalt(isnan(BBaltRaw))=1;
            diffBBalt=diff(maskBBalt);
            
            newIndsMask=1:length(BBaltRaw);
            
            endInds=find(diffBBalt==-1);
            startInds=find(diffBBalt==1);
            
            startInds=startInds+1;
            
            if endInds(1)<startInds(1)
                startInds=[1 startInds];
            end
            if length(endInds)~=length(startInds)
                endInds=[endInds length(maskBBalt)];
            end
            
            % Go through each nan stretch
            for ll=1:length(startInds)
                nanLength=endInds(ll)-startInds(ll);
                % Short stretches
                if nanLength<transLength*2 & startInds(ll)~=1 & endInds(ll)~=length(maskBBalt)
                    BBaltInterp(startInds(ll):endInds(ll))=interp1([startInds(ll)-1,endInds(ll)+1],...
                        [BBaltRaw(startInds(ll)-1),BBaltRaw(endInds(ll)+1)],startInds(ll):endInds(ll));
                else
                    % Tail
                    if startInds(ll~=1)
                        startTail=startInds(ll);
                        endTail=min([startInds(ll)+transLength,length(maskBBalt)]);
                        modT=layerAltsTemp(startTail:endTail);
                        modTtail=layerAltsTemp(endTail)-(layerAltsTemp(endTail)-layerAltsTemp(startTail));
                        int1=interp1([startTail-1,endTail+1],...
                            [BBaltRaw(startTail-1),modTtail],startTail:endTail);
                        addInt=int1-int1(end);
                        BBaltZero(startTail:endTail)=modT+addInt;
                    end
                    % Head
                    if endInds(ll)~=length(maskBBalt)
                        startHead=max([endInds(ll)-transLength,1]);
                        endHead=endInds(ll);
                        modT=layerAltsTemp(startHead:endHead);
                        modThead=layerAltsTemp(startHead)-(layerAltsTemp(startHead)-layerAltsTemp(endHead));
                        int1=interp1([startHead-1,endHead+1],...
                            [modThead,BBaltRaw(endHead+1)],startHead:endHead);
                        addInt=int1-int1(1);
                        BBaltZero(startHead:endHead)=modT+addInt;
                    end
                end
            end
            
            BBaltZero(isnan(BBaltZero))=layerAltsTemp(isnan(BBaltZero));
            BBaltZero(~isnan(BBaltRaw))=nan;
            BBaltZero(~isnan(BBaltInterp))=nan;
        else
            BBaltZero=layerAltsAdj(kk,timeInds);
        end
        
        BBaltAll{end+1}=BBaltRaw;
        timeIall{end+1}=timeInds;
        BBaltInterpAll{end+1}=BBaltInterp;
    else
        BBaltZero=layerAltsAdj(kk,timeInds);
    end
    BBaltZeroAll{end+1}=BBaltZero;
    timeIZeroAll{end+1}=timeInds;
end

%% Put things back together
BBfinished=nan(size(data.LDR));

% Zero degree alt
for ii=1:size(layerAlts,1)
    rows1=layerInds(ii,:);
    cols1=1:length(data.time);
    
    indsNonNan=cat(2,rows1',cols1');
    indsNonNan(any(isnan(indsNonNan), 2), :) = [];
    
    linInds1=sub2ind(size(data.asl),indsNonNan(:,1),indsNonNan(:,2));
    BBfinished(linInds1)=0;
end

% BB altitude
for ii=1:size(BBaltAll,2)
    BB1=BBaltAll{ii};
    timeI1=timeIall{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(data.asl(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=1;
    end
end

% BB altitude interpolated
for ii=1:size(BBaltInterpAll,2)
    BB1=BBaltInterpAll{ii};
    timeI1=timeIall{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(data.asl(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=2;
    end
end

% BB altitude zero
for ii=1:size(BBaltZeroAll,2)
    BB1=BBaltZeroAll{ii};
    timeI1=timeIZeroAll{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(data.asl(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=3;
    end
end

end