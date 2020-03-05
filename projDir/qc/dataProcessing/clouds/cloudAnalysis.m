% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,24,0,0,0);
% endTime=datetime(2018,2,24,23,0,0);

startTime=datetime(2019,8,7,0,0,0);
endTime=datetime(2019,8,7,23,0,0);

plotFields=0;
plotWholeFlight=1;

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
upInd=find(data.elevation>0);

tempVec=1:1:size(data.TEMP,1);
tempMat=repmat(tempVec,size(data.TEMP,2),1)';

tempMat(data.TEMP>0)=nan;
tempMat(isnan(data.TEMP))=nan;

[maxs,rowInds] = nanmax(tempMat);
colInds=1:1:size(data.TEMP,2);

linearInd = sub2ind(size(data.TEMP), rowInds, colInds);

[mins,rowIndsUp] = nanmin(tempMat);
colIndsUp=1:1:size(data.TEMP,2);

linearIndUp = sub2ind(size(data.TEMP), rowIndsUp, colIndsUp);

linearInd(upInd)=linearIndUp(upInd);
rowInds(upInd)=rowIndsUp(upInd);

checkTemps=data.TEMP(linearInd);
outInds=zeros(size(checkTemps));
outInds(checkTemps>0 | checkTemps<-2)=1;
outInds(isnan(checkTemps))=1;
linearInd(outInds==1)=[];
zTempsNeg=data.asl(linearInd);
[zSubR zSubC]=ind2sub(size(data.asl),linearInd);

zTnegFull=nan(size(data.time));
zTnegFull(zSubC)=zTempsNeg;

timeMat=repmat(data.time,size(data.TEMP,1),1);

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
    plot(timeMat(linearInd),zTempsNeg./1000,'c','linewidth',2);
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
%% Find bright band

disp('Searching melting layer ...');
maxLevel=nan(size(data.time));

data.LDR(data.FLAG>1)=nan;
BB=data.LDR;
BB(BB<-16 | BB>-7)=nan;

BB=smoothdata(BB,1,'movmedian',5);

% Find altitude of bright band
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

BBaltRaw=data.asl(linearIndBB);
BBalt=movmedian(BBaltRaw,300,'omitnan');

BBaltFull=nan(size(data.time));
BBaltFull(find(~isnan(maxLevel)))=BBalt;

BBloc=abs(BBaltRaw-BBalt);

BBlocFull=nan(size(data.time));
BBlocFull(find(~isnan(maxLevel)))=BBloc;

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
            %if minVal>-14
            BB((BBvalFull(ii)-BBcol)>(BBvalFull(ii)-minVal)/1.2,ii)=nan;
            %end
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

%% Plot
if plotFields
    fig2=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,800]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    
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
    
    %% Re-sample altitude
    BBdata=sum(BB,1,'omitnan');
    BBnan=find(BBdata==0);
    
    BBaltFull(BBnan)=nan;
    
    newInds=1:100:length(data.time);
    newDBZ=data.DBZ(:,newInds);
    newASL=data.asl(:,newInds);
    newTime=data.time(newInds);
    
    fig3=figure('DefaultAxesFontSize',11,'position',[100,1300,2500,600]);
    
    hold on;
    sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(data.time,BBaltFull./1000,10,'c','filled');
    plot(data.time,zTnegFull./1000,'-k','linewidth',1.5);
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'wholeFlight_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
end