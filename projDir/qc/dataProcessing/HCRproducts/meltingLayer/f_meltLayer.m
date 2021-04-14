% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function [BBfinishedOut iceLevAsl offset]= f_meltLayer(data,zeroAdjustMeters)

debugFig=0;

%% Find zero degree altitude

disp('Searching zero degree altitude ...');

oneGate=data.range(2)-data.range(1);
zeroAdjustGates=round(zeroAdjustMeters/oneGate);

tempTemp=data.TEMP;
%tempTemp(1:16,:)=nan;
oddAngles=find(data.elevation>-70 & data.elevation<70);
tempTemp(:,oddAngles)=nan;
signTemp=sign(tempTemp);
zeroDeg=diff(signTemp,1);

clear tempTemp signTemp

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(zeroDeg~=0)=1;

zeroSum=sum(zeroDeg,1);
tempDataY=find(zeroSum~=0);
zeroSum=zeroSum(tempDataY);
zeroDeg=zeroDeg(:,tempDataY);

zeroAlt=nan(size(zeroDeg));
aslTemp=data.asl(:,tempDataY);
zeroAlt(zeroDeg==1)=aslTemp(zeroDeg==1);

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

clear zeroTemp

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
            elseif minDiff~=100000;
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        end
    end
end

clear zeroAlt zeroDeg

%% Add adjusted layers

disp('Adjusting zero degree layers ...');
layerIndsAdj=nan(size(layerInds));
layerAltsAdj=nan(size(layerAlts));

adjustMeters=zeroAdjustGates*oneGate;
elevTemp=data.elevation(tempDataY);

for kk=1:length(tempDataY)
    if elevTemp(kk)<0
        layerIndsAdj(:,kk)=layerInds(:,kk)-zeroAdjustGates;
    else
        layerIndsAdj(:,kk)=layerInds(:,kk)+zeroAdjustGates;
    end
    layerAltsAdj(:,kk)=layerAlts(:,kk)+adjustMeters;
end

%% Remove data that is not suitable
LDRdata=data.LDR(:,tempDataY);
tempFlag=data.FLAG(:,tempDataY);
%data=rmfield(data,'FLAG');
LDRdata(tempFlag>1)=nan;
LDRdata(LDRdata<-16 | LDRdata>-7)=nan;
LDRdata(1:20,:)=nan;

VELdata=data.VEL_CORR(:,tempDataY);
VELdata(tempFlag>1)=nan;
VELdata(1:20,:)=nan;

clear tempFlag

%% Tightened backlobe
% Initiate mask
blMask=zeros(size(aslTemp));

tempDBZ=data.DBZ(:,tempDataY);
tempWIDTH=data.WIDTH(:,tempDataY);
blMask(tempDBZ<-18 & tempWIDTH>1)=1;

clear tempWIDTH tempDBZ

% Remove small areas
%blMask(:,find(elevTemp<0))=1;
blMask=bwareaopen(blMask,10);

% Fill holes
blMask=imfill(blMask,'holes');

% Only when scanning up
blMask(:,find(elevTemp<0))=0;

% Only within right altitude
rightAlt=(data.altitude(tempDataY)-data.TOPO(tempDataY))*2;

altMat=repmat(rightAlt,size(data.asl(:,tempDataY),1),1);
% Lower limit
blMask(data.asl(:,tempDataY)<(altMat-600))=0;
% Upper limit
blMask(data.asl(:,tempDataY)>(altMat+600))=0;

LDRdata(blMask==1)=nan;
VELdata(blMask==1)=nan;

clear blMask

%% Find altitude of bright band

disp('Searching for melting layer ...');

BBaltAll={};
timeIall={};
BBaltInterpAll={};
BBaltZeroAll={};
timeIZeroAll={};

% Transition length should be 600 seconds
resol=seconds(median(diff(data.time)));
transLength=1/resol*600;

BBaltRawAll=cell(1,size(layerAltsAdj,1));

for kk=1:size(layerAltsAdj,1)
    timeInds=find(~isnan(layerAltsAdj(kk,:)) & layerIndsAdj(kk,:)>0);
    layerAltsTemp=layerAltsAdj(kk,timeInds);
    rowInds=layerIndsAdj(kk,timeInds);
    
    % LDR
    % Find maximum LDR level in area around zero degree altitude
    maxLevelLDR=nan(size(rowInds));
    for ii=1:length(rowInds)
        vertColLDR=LDRdata(:,timeInds(ii));
        vertColLDR(rowInds(ii)+18:end)=nan;
        vertColLDR(1:max([rowInds(ii)-18,1]))=nan;
        
        % Check if all nan
        if min(isnan(vertColLDR))==0
            %vertColData=vertColLDR(~isnan(vertColLDR));
            maxLevelLDR(ii)=min(find(vertColLDR==nanmax(vertColLDR)));
        end
    end
    
    ASLlayer=aslTemp(:,timeInds);
    
    if min(isnan(maxLevelLDR))==0
        
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
        
        % Remove data where distance is more than 50 m
        LDRaltRaw(LDRloc>50)=nan;
        
        % Adjust zero degree layer
        zeroDist=LDRaltRaw-layerAltsTemp;
        
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
                        
            distInds=round(noNanDist./oneGate);
            rowInds(elevTemp(timeInds)>0)=rowInds(elevTemp(timeInds)>0)-distInds(elevTemp(timeInds)>0);
            rowInds(elevTemp(timeInds)<=0)=rowInds(elevTemp(timeInds)<=0)+distInds(elevTemp(timeInds)<=0);
        end
    end
    
    % VEL
    
    % Find vel melting layer level
    maxLevelVEL=nan(size(rowInds));
    velMasked=VELdata;
    velMasked(:,elevTemp<0)=-velMasked(:,elevTemp<0);
    
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
    
    % Find altitude of step
    stepInLin=find(velSteps==1);
    [stepInR,stepInC]=ind2sub(size(velSmooth),stepInLin);
    
    maxLevelVEL=nan(1,length(timeInds));
    maxLevelVEL(stepInC)=stepInR;
    
    % Remove data that doesn't cut it
    if min(isnan(maxLevelLDR))==0 | min(isnan(maxLevelVEL))==0
                
        % VEL
        colIndsVEL=1:1:size(ASLlayer,2);
        linearIndVEL = sub2ind(size(ASLlayer), maxLevelVEL, colIndsVEL);
        linearIndVEL(isnan(maxLevelVEL))=[];
        
        % Raw altitude of melting layer
        VELaltRawIn=ASLlayer(linearIndVEL);
        VELaltRaw=nan(size(rowInds));
        VELaltRaw(find(~isnan(maxLevelVEL)))=VELaltRawIn;
        
        % Difference between steps
        diffS1=nan(1,size(velSteps,2));
        udS1=nan(2,size(velSteps,2));
        
        for ii=1:size(velSteps,2)
            uS1=unique(S1(:,ii),'stable');
            uS1(isnan(uS1))=[];
            
            if elevTemp(ii)>0
                uS1=flipud(uS1);
            end
            
            if length(uS1)==2
                diffS1(ii)=uS1(2)-uS1(1);
                udS1(:,ii)=uS1;
            end
            
        end
                
        if debugFig & ~isempty(maxLevelVEL)
            fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);
            colmap=jet;
            colormap(flipud(colmap));
            subplot(4,1,1)
            hold on
            surf(velSmooth,'edgecolor','none');
            view(2)
            %colorbar
            caxis([-5 5])
            plot(maxLevelVEL,'-b')
            
            subplot(4,1,2)
            hold on
            plot(diffS1)
            plot(movmean(diffS1,50,'omitnan'),'linewidth',2)
            %ylim([0 3]);
            grid on
        end
        
        VELaltRaw(movmean(udS1(1,:),50,'omitnan')>0)=nan;
        VELaltRaw(movmean(udS1(2,:),50,'omitnan')>0)=nan;
        VELaltRaw(movmean(diffS1,50,'omitnan')>-0.7)=nan;
        
        % Compare with zero degree layer
        VELzeroDiff=VELaltRaw-layerAltsTemp;
        VELaltRaw(abs(VELzeroDiff)>200)=nan;
        
        % Mean altitude of meltig layer
        VELalt=movmedian(VELaltRaw,1000,'omitnan');
        VELalt(isnan(VELaltRaw))=nan;
        
        % Distance between raw and mean altitude
        VELloc=abs(VELaltRaw-VELalt);
        % Remove data where distance is more than 100 m
        VELaltRaw(VELloc>100)=nan;
        
        % Standard deviation
        VELaltS=movstd(VELaltRaw,300,'omitnan');
        VELaltS(isnan(VELaltRaw))=nan;
        % Remove data with too much std
        VELaltRaw(VELaltS>35)=nan;
        
        if debugFig & ~isempty(maxLevelVEL)
            subplot(4,1,3:4)
            hold on
            plot(VELzeroDiff);
            plot(VELloc);
            plot(VELaltRaw,'linewidth',2)
            yyaxis right
            plot(VELaltS)
            legend('Difference to zero deg','Difference from mean','Good data','Std',...
                'location','northoutside','orientation','horizontal');
            grid on
        end
        
        % Combine max levels
        BBaltRaw=VELaltRaw;
        if exist('LDRaltRaw')
            BBaltRaw(~isnan(LDRaltRaw))=LDRaltRaw(~isnan(LDRaltRaw));
        end
        
        clear LDRaltRaw
        
        % Remove outliers
        BBaltRawTest=movmedian(BBaltRaw,30,'includenan');
        BBaltRawTest=fillmissing(BBaltRawTest,'nearest');
        
        BBrawMinusTest=abs(BBaltRawTest-BBaltRaw);
        BBaltRaw(BBrawMinusTest>50)=nan;
        
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
        
        % Remove data that is more than 100 m above model zero degree
        zeroAltTest=layerAlts(kk,timeInds);        
        BBaltCompare=zeroAltTest-BBaltRaw;
        BBaltRaw(BBaltCompare<-100)=nan;
        
        BBaltRawAll{kk}=BBaltRaw;
    end
end

% Calculate mean offset for whole flight
disp('Calculating offset ...');
zeroDistCollect=[];
offset=[];
for kk=1:size(layerAlts,1)
    timeInds=find(~isnan(layerAltsAdj(kk,:)) & layerIndsAdj(kk,:)>0);
    BBaltRaw=BBaltRawAll{kk};
    if ~isempty(BBaltRaw)
        zeroDistAll=BBaltRaw-layerAlts(kk,timeInds);
        zeroDistAll(isnan(zeroDistAll))=[];
        zeroDistCollect=[zeroDistCollect zeroDistAll];
    end
end

if ~isnan(zeroDistCollect)
    offset=median(zeroDistCollect);
end

% Adjust zero degree layer
for kk=1:size(layerAltsAdj,1)
    if ~isempty(offset) & ~isempty(BBaltRawAll{kk})
        BBaltRaw=BBaltRawAll{kk};
        timeInds=find(~isnan(layerAltsAdj(kk,:)) & layerIndsAdj(kk,:)>0);
        zeroDistF=BBaltRaw-layerAlts(kk,timeInds);
        
        % Remove small data stretches
        zeroMaskF=zeros(size(zeroDistF));
        zeroMaskF(~isnan(zeroDistF))=1;
        
        zeroMaskF=bwareaopen(zeroMaskF,5);
        zeroDistF(zeroMaskF==0)=nan;
        
        % Adjust end points so that they have the median distance
        firstInd=min(find(~isnan(zeroDistF)));
        adjPixNumStart=round(max([1 firstInd-1000]));
        zeroDistF(1:adjPixNumStart)=offset;
        lastInd=max(find(~isnan(zeroDistF)));
        adjPixNumEnd=round(min([length(zeroDistF) lastInd+1000]));
        zeroDistF(adjPixNumEnd:end)=offset;
        
        zeroDistSmoothF=movmean(zeroDistF,3000,'omitnan');
        zeroDistSmoothF(isnan(zeroDistF))=nan;
        
        noNanDistF=fillmissing(zeroDistSmoothF,'linear','EndValues','nearest');
        layerAltsAdj(kk,timeInds)=layerAlts(kk,timeInds)+noNanDistF;
    end
end

%% Interpolation
disp('Interpolating ...');

for kk=1:size(layerAltsAdj,1)
    BBaltRaw=BBaltRawAll{kk};
    timeInds=find(~isnan(layerAltsAdj(kk,:)) & layerIndsAdj(kk,:)>0);
    layerAltsTempF=layerAltsAdj(kk,timeInds);
        
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
        
        if isempty(startInds) | endInds(1)<startInds(1)
            startInds=[1 startInds];
        end
        if length(endInds)~=length(startInds)
            endInds=[endInds length(maskBBalt)];
        end
        
        % Go through each nan stretch
        for ll=1:length(startInds)
            nanLength=endInds(ll)-startInds(ll);
            % Short stretches
            if nanLength<transLength*2 & startInds(ll)~=1 & endInds(ll)~=length(maskBBalt) ...
                    & abs(BBaltRaw(startInds(ll)-1)-BBaltRaw(endInds(ll)+1))<500
                BBaltInterp(startInds(ll):endInds(ll))=interp1([startInds(ll)-1,endInds(ll)+1],...
                    [BBaltRaw(startInds(ll)-1),BBaltRaw(endInds(ll)+1)],startInds(ll):endInds(ll));
            else
                % Tail
                if startInds(ll~=1)
                    startTail=startInds(ll);
                    endTail=min([startInds(ll)+transLength,length(maskBBalt)]);
                    modT=layerAltsTempF(startTail:endTail);
                    modTtail=layerAltsTempF(endTail)-(layerAltsTempF(endTail)-layerAltsTempF(startTail));
                    int1=interp1([startTail-1,endTail+1],...
                        [BBaltRaw(startTail-1),modTtail],startTail:endTail);
                    addInt=int1-int1(end);
                    BBaltZero(startTail:endTail)=modT+addInt;
                end
                % Head
                if endInds(ll)~=length(maskBBalt)
                    startHead=max([endInds(ll)-transLength,1]);
                    endHead=endInds(ll);
                    modT=layerAltsTempF(startHead:endHead);
                    modThead=layerAltsTempF(startHead)-(layerAltsTempF(startHead)-layerAltsTempF(endHead));
                    int1=interp1([startHead-1,endHead+1],...
                        [modThead,BBaltRaw(endHead+1)],startHead:endHead);
                    addInt=int1-int1(1);
                    BBaltZero(startHead:endHead)=modT+addInt;
                end
            end
        end
        
        BBaltZero(isnan(BBaltZero))=layerAltsTempF(isnan(BBaltZero));
        BBaltZero(~isnan(BBaltRaw))=nan;
        BBaltZero(~isnan(BBaltInterp))=nan;
    else
        BBaltZero=layerAltsAdj(kk,timeInds);
    end
    
    BBaltAll{end+1}=BBaltRaw;
    timeIall{end+1}=timeInds;
    BBaltInterpAll{end+1}=BBaltInterp;
    BBaltZeroAll{end+1}=BBaltZero;
    timeIZeroAll{end+1}=timeInds;
end

clear LDRdata VELdata

%% Put things back together
BBfinished=nan(size(aslTemp));

disp('Creating melting layer output ...');

% Zero degree alt (0)
for ii=1:size(layerAlts,1)
    rows1=layerInds(ii,:);
    cols1=1:length(tempDataY);
    
    indsNonNan=cat(2,rows1',cols1');
    indsNonNan(any(isnan(indsNonNan), 2), :) = [];
    
    linInds1=sub2ind(size(aslTemp),indsNonNan(:,1),indsNonNan(:,2));
    BBfinished(linInds1)=0;
end

% BB altitude (1)
for ii=1:size(BBaltAll,2)
    BB1=BBaltAll{ii};
    timeI1=timeIall{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(aslTemp(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=1;
    end
end

% BB altitude interpolated (2)
for ii=1:size(BBaltInterpAll,2)
    BB1=BBaltInterpAll{ii};
    timeI1=timeIall{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(aslTemp(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=2;
    end
end

% BB altitude inferred from zero (3)
for ii=1:size(BBaltZeroAll,2)
    BB1=BBaltZeroAll{ii};
    timeI1=timeIZeroAll{ii};
    
    timeI1(isnan(BB1))=[];
    BB1(isnan(BB1))=[];
    
    for jj=1:length(BB1)
        [min1 ind1]=(min(abs(aslTemp(:,timeI1(jj))-BB1(jj))));
        BBfinished(ind1,timeI1(jj))=3;
    end
end

BBfinishedOrigInds=nan(size(data.DBZ));
BBfinishedOrigInds(:,tempDataY)=BBfinished;

clear BBfinished;
BBfinishedOrigInds(data.asl<0)=nan;
BBfinishedOrigInds(1,:)=nan;

%% Change assignments and create icing level
% Below icing level=10
% Above icing level=20
% Melting or 0 deg layer at or below icing level=11 (0 deg), 12 (detected),
% 13 (interpolated), or 14 (estimated)
% Melting or 0 deg layer above icing level=21, 22, 23, or 24

disp('Finding icing level ...');
BBfinishedOut=nan(size(BBfinishedOrigInds));
iceLev=nan(1,size(BBfinishedOut,2));

for ii=1:size(BBfinishedOrigInds,2)
    BBcol=BBfinishedOrigInds(:,ii);
    altCol=data.asl(:,ii);
    oneAlts=altCol(find(BBcol>0));
    if ~isempty(oneAlts)
        iceLevInd=find(altCol==min(oneAlts));
        iceLev(ii)=data.asl(iceLevInd,ii);
        if data.elevation(ii)<0 % Pointing down
            BBfinishedOut(1:iceLevInd-1,ii)=20;
            BBfinishedOut(iceLevInd:end,ii)=10;
        else % Pointing up
            BBfinishedOut(1:iceLevInd,ii)=10;
            BBfinishedOut(iceLevInd+1:end,ii)=20;
        end
    end
end

% Count number of zero deg and melt layers

numZeroDeg=sum(BBfinishedOrigInds==0,1);
numOther=sum(BBfinishedOrigInds>0,1);

diffNum=numZeroDeg-numOther;

% Remove data with odd angles
iceLev(oddAngles)=nan;

% Take care of areas with big jumps
%[change1 S1]=ischange(iceLev,'linear','Threshold',1000000);

% Remove outliers
iceLevMed=movmedian(iceLev,10,'includenan');
iceLevMed(isnan(iceLev))=nan;
iceLev(abs(iceLevMed-iceLev)>100)=nan;

iceLevDiff=iceLev;

iceLevDiff=fillmissing(iceLevDiff,'next');
diff1=diff(iceLevDiff);
change1=zeros(size(iceLev));
change1(abs(diff1)>70)=1;

changeInds=find(change1==1);
changeStart=[1 changeInds];
changeEnd=[changeInds+1 length(iceLev)];

iceLevTest=iceLev;

iceShort1=iceLev(changeStart(1):changeEnd(1));
iceShort1(isnan(iceShort1))=[];
currentLev=median(iceShort1(max([1, length(iceShort1)-200]):end),'omitnan');

for ii=1:length(changeStart)
    iceShort=iceLev(changeStart(ii):changeEnd(ii));
    iceShort(isnan(iceShort))=[];
    changeLength=length(iceShort);
    newLev=median(iceShort(1:min([200,length(iceShort)])),'omitnan');
    
    moreNum=diffNum(changeStart(ii):changeEnd(ii));
    sumMoreNum=length(find(moreNum>0));
    
    if (changeLength<1000 & abs(newLev-currentLev)>100) | (changeLength<5000 & sumMoreNum>50 & abs(newLev-currentLev)>100)
        iceLevTest(changeStart(ii):changeEnd(ii))=nan;
    else
        currentLev=median(iceShort(max([1, length(iceShort)-200]):end),'omitnan');
    end
end

iceLevTest=fillmissing(iceLevTest,'linear','endValues','nearest');
iceLevDiff=abs(iceLev-iceLevTest);
iceLev(iceLevDiff>100)=nan;

missInds=find(isnan(iceLev));

smoothIce=iceLev;
smoothIce=fillmissing(smoothIce,'linear','endValues','nearest');

% Fill in missing
for ii=1:length(missInds)
    BBfinishedOut(find(data.asl(:,missInds(ii))<smoothIce(missInds(ii))),missInds(ii))=10;
    BBfinishedOut(find(data.asl(:,missInds(ii))>=smoothIce(missInds(ii))),missInds(ii))=20;
end

BBfinishedOut(:,oddAngles)=nan;
BBfinishedOut(data.asl<=0)=nan;

iceLevAslAll=data.asl;
iceLevAslAll(BBfinishedOut~=20)=nan;
[minIce minInd]=min(iceLevAslAll,[],1);
linInd=sub2ind(size(data.asl),minInd,1:size(data.asl,2));
iceLevAsl=data.asl(linInd);

% Remove nan
iceLevAsl(find(all(isnan(iceLevAslAll),1)))=nan;

% Remove data with no 10s
iceLevAslAll=data.asl;
iceLevAslAll(BBfinishedOut~=10)=nan;
iceLevAsl(find(all(isnan(iceLevAslAll),1)))=nan;

% Remove data close to the ground
groundInds=find(iceLevAsl<30 & isnan(iceLev));
iceLevAsl(groundInds)=nan;
for ii=1:length(groundInds)
    BBfinishedOut(find(BBfinishedOut(:,groundInds(ii))==10),groundInds(ii))=20;
end

% Zero deg
BBfinishedOut(BBfinishedOrigInds==0 & BBfinishedOut==10)=11;
BBfinishedOut(BBfinishedOrigInds==0 & BBfinishedOut==20)=21;

% detected
BBfinishedOut(BBfinishedOrigInds==1 & BBfinishedOut==10)=12;
BBfinishedOut(BBfinishedOrigInds==1 & BBfinishedOut==20)=22;

% interpolated
BBfinishedOut(BBfinishedOrigInds==2 & BBfinishedOut==10)=13;
BBfinishedOut(BBfinishedOrigInds==2 & BBfinishedOut==20)=23;

% estimated
BBfinishedOut(BBfinishedOrigInds==3 & BBfinishedOut==10)=14;
BBfinishedOut(BBfinishedOrigInds==3 & BBfinishedOut==20)=24;

if isempty(offset)
    offset=nan;
end
end