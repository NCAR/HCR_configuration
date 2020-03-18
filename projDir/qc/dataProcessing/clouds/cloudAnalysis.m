% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2018,2,4,12,10,0);
endTime=datetime(2018,2,5,12,35,0);

% startTime=datetime(2019,8,7,0,0,0);
% endTime=datetime(2019,8,7,23,0,0);

plotFields=0;
plotWholeFlight=1;

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

% Expected bright band altitude. Determines plot zoom.
if strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 15];
elseif strcmp(project,'socrates')
    expBBalt=2;
    ylimits=[-0.2 9];
elseif strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/clouds/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data.DBZ=[];
data.LDR=[];
data.TEMP=[];
data.WIDTH=[];
data.FLAG=[];
data.TOPO=[];
data.pitch=[];
data.roll=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Find zero degree altitude

disp('Searching 0 deg altitude ...');

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

% Remove contiguous layers that are too small
if length(data.time)>9000
    zeroCut=3000;
    CC = bwconncomp(zeroTemp);
    
    for ii=1:CC.NumObjects
        area=CC.PixelIdxList{ii};
        if length(area)<=zeroCut
            zeroTemp(area)=nan;
        end
    end
end

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

%% Remove data that is not suitable
BB=data.LDR;
BB(data.FLAG>1)=nan;
BB(BB<-16 | BB>-7)=nan;
BB(data.range<150)=nan;

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

BB(blMask==1)=nan;

%% Find altitude of bright band

BBaltAll={};
timeIall={};
BBaltInterpAll={};
BBaltZeroAll={};
timeIZeroAll={};

transLength=6000;

for kk=1:size(layerAlts,1)
    timeInds=find(~isnan(layerAlts(kk,:)));
    layerAltsTemp=layerAlts(kk,timeInds);
    
    rowInds=layerInds(kk,timeInds);
    maxLevel=nan(size(rowInds));
    for ii=1:length(rowInds)
        vertCol=BB(:,timeInds(ii));
        if data.elevation(timeInds(ii))<90 % Down pointing
            vertCol(rowInds(ii)+30:end)=nan;
            vertCol(1:max([rowInds(ii)-10,1]))=nan;
        else % Up pointing
            vertCol(rowInds(ii)+10:end)=nan;
            vertCol(1:max([rowInds(ii)-30,1]))=nan;
        end
        % Check if all nan
        if min(isnan(vertCol))==0
            vertColData=vertCol(~isnan(vertCol));
            maxLevel(ii)=min(find(vertCol==nanmax(vertCol)));
        end
    end
    
    if min(isnan(maxLevel))==0
        
        BBlayer=BB(:,timeInds);
        ASLlayer=data.asl(:,timeInds);
        
        colIndsBB=1:1:size(BBlayer,2);
        linearIndBB = sub2ind(size(BBlayer), maxLevel, colIndsBB);
        linearIndBB(isnan(maxLevel))=[];
        
        % Raw altitude of melting layer
        BBaltRawIn=ASLlayer(linearIndBB);
        BBaltRaw=nan(size(rowInds));
        BBaltRaw(find(~isnan(maxLevel)))=BBaltRawIn;
        % Mean altitude of meltig layer
        BBalt=movmedian(BBaltRaw,300,'omitnan');
        BBalt(isnan(BBaltRaw))=nan;
        % Standard deviation
        BBaltS=movstd(BBaltRaw,300,'omitnan');
        BBaltS(isnan(BBaltRaw))=nan;
        % Remove data with too much std
        BBaltRaw(BBaltS>100)=nan;
        % Distance between raw and mean altitude
        BBloc=abs(BBaltRaw-BBalt);
        
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
        
        zeroCut=100;
        CC = bwconncomp(BBmask);
        
        for ii=1:CC.NumObjects
            area=CC.PixelIdxList{ii};
            if length(area)<=zeroCut
                BBaltRaw(area)=nan;
            end
        end
        
        % Remove data where distance is more than 200 m
        BBaltRaw(BBloc>50)=nan;
        
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
                        endTailIn=min([startInds(ll)+transLength,length(maskBBalt)]);
                        modT=layerAltsTemp(startTail:endTailIn);
                        [~,modTind]=min(abs(BBaltRaw(startTail-1)-modT));
                        indsTemp=newIndsMask(startTail:endTailIn);
                        endTail=indsTemp(modTind);
                        
                        modTtail=layerAltsTemp(endTail);
                        BBaltZero(startTail:endTail)=interp1([startTail-1,endTail+1],...
                            [BBaltRaw(startTail-1),modTtail],startTail:endTail);
                    end
                    % Head
                    if endInds(ll)~=length(maskBBalt)
                        startHeadIn=max([endInds(ll)-transLength,1]);
                        endHead=endInds(ll);
                        modT=layerAltsTemp(startHeadIn:endHead);
                        [~,modTind]=min(abs(BBaltRaw(endHead+1)-modT));
                        indsTemp=newIndsMask(startHeadIn:endHead);
                        startHead=indsTemp(modTind);
                        
                        modThead=layerAltsTemp(startHead);
                        BBaltZero(startHead:endHead)=interp1([startHead-1,endHead+1],...
                            [modThead,BBaltRaw(endHead+1)],startHead:endHead);
                    end
                end
            end
            
            BBaltZero(isnan(BBaltZero))=layerAltsTemp(isnan(BBaltZero));
            BBaltZero(~isnan(BBaltRaw))=nan;
            BBaltZero(~isnan(BBaltInterp))=nan;
        else
            BBaltZero=layerAlts(kk,timeInds);
        end
        
        if min(min(isnan(BBlayer)))==0
            BBaltAll{end+1}=BBaltRaw;
            timeIall{end+1}=timeInds;
            BBaltInterpAll{end+1}=BBaltInterp;
        end
    else
        BBaltZero=layerAlts(kk,timeInds);
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

zeroInds=find(BBfinished==0);
oneInds=find(BBfinished==1);
twoInds=find(BBfinished==2);
threeInds=find(BBfinished==3);
%% Plot

close all

if plotFields
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    
    sub1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'k','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
    ax2=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    sub2=colMapDBZ(sub2);
    scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'c','filled');
    scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'b','filled');
    scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'refl_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

%% Plot
if plotFields
    fig2=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,800]);
    
    colormap jet
    %%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    sub2=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
    view(2);
    caxis([-25 5]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    if max(size(BBaltAll))>0
        ax2=subplot(2,1,2);
        hold on;
        sub2=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
        view(2);
        caxis([-25 5]);
        colorbar
        scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'r','filled');
        scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'m','filled');
        scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
        ax = gca;
        ax.SortMethod = 'childorder';
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('LDR and melting layer')
        grid on
    end
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'bb_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

%% DBZ and line plot
if plotWholeFlight
    close all
    
    % Resample for plotting
    newInds=1:100:length(data.time);
    newDBZ=data.DBZ(:,newInds);
    newLDR=data.LDR(:,newInds);
    newASL=data.asl(:,newInds);
    newTEMP=data.TEMP(:,newInds);
    newTime=data.time(newInds);
    
    fig3=figure('DefaultAxesFontSize',11,'position',[100,1300,2500,1000]);
    
    subplot(2,1,1)
    hold on;
    sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    colIndsAll=1:length(data.time);
    for ii=1:size(layerAlts,1)
        rowInds=layerInds(ii,:);
        colInds=colIndsAll;
        colInds(isnan(rowInds))=[];
        rowInds(isnan(rowInds))=[];
        
        linPlot=sub2ind(size(data.asl),rowInds,colInds);
        
        rowAlts=layerAlts(ii,:);
        rowAlts(isnan(rowAlts))=[];
        scatter(timeMat(linPlot),rowAlts./1000,10,'filled');
    end
    scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'c','filled');
    scatter(timeMat(twoInds),data.asl(twoInds)./1000,10,'b','filled');
    scatter(timeMat(threeInds),data.asl(threeInds)./1000,10,'g','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
        
    subplot(2,1,2)
    hold on;
    sub2=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
    view(2);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'wholeFlight_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
end