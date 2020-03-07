% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2018,2,7,18,0,0);
endTime=datetime(2018,2,8,12,0,0);

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
zeroCut=3000;
CC = bwconncomp(zeroTemp);

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=zeroCut
        zeroTemp(area)=nan;
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
BB(data.range<500)=nan;

% Smooth data along time axis
%BB=smoothdata(BB,1,'movmedian',5);

% Find altitude of bright band
maxLevel=nan(size(data.time));
for ii=1:length(data.time)
    vertCol=BB(:,ii);
    if data.elevation(ii)<90 % Down pointing
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

colIndsBB=1:1:size(data.TEMP,2);
linearIndBB = sub2ind(size(data.TEMP), maxLevel, colIndsBB);
linearIndBB(isnan(maxLevel))=[];

% Raw altitude of melting layer
BBaltRaw=data.asl(linearIndBB);
% Mean altitude of meltig layer
BBalt=movmedian(BBaltRaw,300,'omitnan');
% Distance between raw ans mean altitude
BBloc=abs(BBaltRaw-BBalt);

% Resample to full resolution
BBaltFull=nan(size(data.time));
BBaltFull(find(~isnan(maxLevel)))=BBalt;
BBlocFull=nan(size(data.time));
BBlocFull(find(~isnan(maxLevel)))=BBloc;

% Remove data where distance is more than 200 m
BB(:,isnan(BBlocFull))=nan;
BB(:,BBlocFull>200)=nan;

%% Remove data that is too far off bright band altitude

BBvalMax=BB(linearIndBB);
BBvalFull=nan(size(data.time));
BBvalFull(find(~isnan(maxLevel)))=BBvalMax;

for ii=1:length(data.time)
    BBcol=BB(:,ii);
    dataInd=find(~isnan(BBcol));
    if length(dataInd)>0
        diffs=abs(data.asl(dataInd,ii)-BBaltFull(ii));
        tooFar=find(diffs>300);
        BB(dataInd(tooFar),ii)=nan;
        % remove data in low gradient areas
        BBcol=BB(:,ii);
        dataInd=find(~isnan(BBcol));
        if length(dataInd)>0
            minVal=min(BBcol(dataInd));
            BB((BBvalFull(ii)-BBcol)>(BBvalFull(ii)-minVal)/1.2,ii)=nan;
        end
    end
end

% Remove contiguous areas that are too small
specCut=100;

bbTemp=BB;

maskTemp=zeros(size(bbTemp));
maskTemp(~isnan(bbTemp))=1;

CC = bwconncomp(maskTemp,4);

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=specCut
        BB(area)=nan;
    end
end

%% Re-sample altitude (we keep only data where BB is nview(2);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);on nan)
BBdata=sum(BB,1,'omitnan');
BBnan=find(BBdata==0);

BBaltFull(BBnan)=nan;

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
    scatter(timeMat(zeroDeg==1),zeroAlt./1000,10,'k','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
    ax2=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
    view(2);
    %plot(timeMat(linearInd),zTempsNeg./1000,'c','linewidth',2);
    colorbar
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'refl_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

%% Plot
if plotFields
    fig2=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,800]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ax1=subplot(2,1,1);
    hold on;
    sub1=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
    view(2);
    scatter(timeMat(linearIndBB),BBalt./1000,10,'m','filled');
    scatter(timeMat(linearIndBB),data.asl(linearIndBB)./1000,10,'r','filled');
    colorbar
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim([expBBalt-3 expBBalt+3]);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    ax2=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,BB,'edgecolor','none');
    view(2);
    colorbar
    ylim([expBBalt-3 expBBalt+3]);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Bright band')
    grid on
    
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
    %scatter(data.time,BBaltFull./1000,10,'c','filled');
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