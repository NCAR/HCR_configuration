% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,8,7,12,35,0);
endTime=datetime(2019,8,7,12,53,0);

plotFields=1;
plotWholeFlight=0;

project='otrec'; %socrates, aristo, cset
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
data.FLAG=[];
data.pitch=[];

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
BB(data.range<100)=nan;

% Smooth data along time axis
%BB=smoothdata(BB,1,'movmedian',5);

%% Find altitude of bright band
BBall={};
BBaltAll={};
timeIall={};

for kk=1:size(layerAlts,1)
    timeInds=find(~isnan(layerAlts(kk,:)));
    
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
    
    if min(isnan(maxLevel))~=0
        continue
    end
    
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
    BBalt(BBaltS>150)=nan;
    % Distance between raw and mean altitude
    BBloc=abs(BBaltRaw-BBalt);
    
    % Remove data where distance is more than 200 m
    BBlayer(:,isnan(BBloc))=nan;
    BBlayer(:,BBloc>200)=nan;
    
    % Remove data that is too far off bright band altitude
    
    BBvalMax=BBlayer(linearIndBB);
    BBvalFull=nan(size(rowInds));
    BBvalFull(find(~isnan(maxLevel)))=BBvalMax;
    
    for ii=1:length(rowInds)
        BBcol=BBlayer(:,ii);
        dataInd=find(~isnan(BBcol));
        if length(dataInd)>0
            diffs=abs(ASLlayer(dataInd,ii)-BBalt(ii));
            tooFar=find(diffs>300);
            BBlayer(dataInd(tooFar),ii)=nan;
            % remove data in low gradient areas
            BBcol=BBlayer(:,ii);
            dataInd=find(~isnan(BBcol));
            if length(dataInd)>0
                minVal=min(BBcol(dataInd));
                BBlayer((BBvalFull(ii)-BBcol)>(BBvalFull(ii)-minVal)/1.2,ii)=nan;
            end
        end
    end
    
    % Remove contiguous areas that are too small
    specCut=100;
    
    bbTemp=BBlayer;
    
    maskTemp=zeros(size(bbTemp));
    maskTemp(~isnan(bbTemp))=1;
    
    CC = bwconncomp(maskTemp,4);
    
    for ii=1:CC.NumObjects
        area=CC.PixelIdxList{ii};
        if length(area)<=specCut
            BBlayer(area)=nan;
        end
    end
    
    % Re-sample altitude (we keep only data where BB has data);
    
    BBdata=sum(BBlayer,1,'omitnan');
    BBnan=find(BBdata==0);
    
    BBalt(BBnan)=nan;
    
    if min(min(isnan(BBlayer)))==0
        BBall{end+1}=BBlayer;
        BBaltAll{end+1}=BBalt;
        timeIall{end+1}=timeInds;
    end
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

zeroInds=find(BBfinished==0);
oneInds=find(BBfinished==1);
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
    ax2=subplot(2,1,1);
    hold on;
    sub2=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
    view(2);
    caxis([-25 5]);
    colorbar
    scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'r','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    if max(size(BBaltAll))>0
        ax2=subplot(2,1,2);
        hold on;
        sub2=surf(data.time,data.asl./1000,BBlayer,'edgecolor','none');
        view(2);
        caxis([-25 5]);
        colorbar
        ylim(ylimits);
        ylabel('Altitude (km)');
        xlim([data.time(1),data.time(end)]);
        title('Bright band')
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
    scatter(timeMat(oneInds),data.asl(oneInds)./1000,10,'c','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
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
    
    %scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'b','filled');
    
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