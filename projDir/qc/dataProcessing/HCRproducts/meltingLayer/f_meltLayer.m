% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function [BBfinishedOut]= f_meltLayer(data,zeroAdjustMeters)

debugFig=0;

%% Find zero degree altitude

disp('Searching zero degree altitude ...');

oneGate=data.range(2)-data.range(1);
zeroAdjustGates=round(zeroAdjustMeters/oneGate);

tempTemp=data.TEMP;
tempTemp(1:16,:)=nan;
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
        layerIndsAdj(:,kk)=layerInds(:,kk)+zeroAdjustGates;
    else
        layerIndsAdj(:,kk)=layerInds(:,kk)-zeroAdjustGates;
    end
    layerAltsAdj(:,kk)=layerAlts(:,kk)-adjustMeters;
end

%% Remove data that is not suitable
LDRdata=data.LDR(:,tempDataY);
tempFlag=data.FLAG(:,tempDataY);
%data=rmfield(data,'FLAG');
LDRdata(tempFlag>1)=nan;
LDRdata(LDRdata<-16 | LDRdata>-7)=nan;
tempRange=data.range(tempDataY);
LDRdata(tempRange<150)=nan;
LDRdata(1:20,:)=nan;

VELdata=data.VEL_CORR(:,tempDataY);
VELdata(tempFlag>1)=nan;
VELdata(tempRange<150)=nan;
VELdata(1:20,:)=nan;

clear tempFlag

%% Tightened backlobe
% Initiate mask
blMask=zeros(size(aslTemp));

tempDBZ=data.DBZ(:,tempDataY);
tempWIDTH=data.WIDTH(:,tempDataY);
blMask(tempDBZ<-18 & tempWIDTH>1)=1;

clear tempWIDTH tempDBZ

% Only within right altitude
rightAlt=data.altitude(tempDataY)-data.TOPO(tempDataY);

altMat=repmat(rightAlt,size(tempRange,1),1);
% Lower limit
blMask(tempRange<(altMat-100))=0;
% Upper limit
blMask(tempRange>(altMat+2500))=0;

% Only when scanning up
blMask(:,find(elevTemp<0))=0;

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
        disp('Searching LDR ...');
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
    disp('Searching VEL ...');
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
        disp('Removing bad data ...');
        
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
            uS1=unique(S1(:,ii));
            uS1(isnan(uS1))=[];
            
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
            ylim([0 3]);
            grid on
        end
        
        VELaltRaw(movmean(udS1(1,:),50,'omitnan')>0)=nan;
        VELaltRaw(movmean(udS1(2,:),50,'omitnan')>0)=nan;
        VELaltRaw(movmean(diffS1,50,'omitnan')<0.7)=nan;
        
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
                
        % Interpolate between good values
        BBaltInterp=nan(size(BBaltRaw));
        BBaltZero=nan(size(BBaltRaw));
        
        if min(isnan(BBaltRaw))==0
            
            % Adjust zero degree layer
            zeroDistF=BBaltRaw-layerAltsAdj(kk,timeInds);
            
            % Remove small data stretches
            zeroMaskF=zeros(size(zeroDistF));
            zeroMaskF(~isnan(zeroDistF))=1;
            
            zeroMaskF=bwareaopen(zeroMaskF,5);
            zeroDistF(zeroMaskF==0)=nan;
            
            zeroDistSmoothF=movmean(zeroDistF,3000,'omitnan');
            zeroDistSmoothF(isnan(zeroDistF))=nan;
            
            noNanDistF=fillmissing(zeroDistSmoothF,'linear','EndValues','nearest');
            if sum(~isnan(noNanDistF))>0
                layerAltsTempF=layerAltsTemp+noNanDistF;
            end
            
            layerAltsAdj(kk,timeInds)=layerAltsTempF;
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
                if nanLength<transLength*2 & startInds(ll)~=1 & endInds(ll)~=length(maskBBalt)
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
    else
        BBaltZero=layerAltsAdj(kk,timeInds);
    end
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

BBfinishedOut=nan(size(data.DBZ));
BBfinishedOut(:,tempDataY)=BBfinished;
end