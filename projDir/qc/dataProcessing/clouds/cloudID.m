% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,10,2,14,50,0);
endTime=datetime(2019,10,2,14,54,0);

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

figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudID/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data.DBZ=[];
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

zeroInds=find(~isnan(zeroAlt));

%% Mask non cloud

refl=data.DBZ;
refl(data.FLAG>1)=nan;

%% Add extinct back in

flagTemp=data.FLAG;
flagTemp(data.FLAG==8)=7;
surfMask=ones(1,size(flagTemp,2));
surfMask(find(any(data.FLAG==7 | data.FLAG==8,1)))=0;
surfDiff=diff(surfMask);

% Add surface echo back in
holeStart=find(surfDiff==1);
holeStart=holeStart+1;
holeEnd=find(surfDiff==-1);

if holeStart>holeEnd
    holeStart=[1 holeStart];
end
if length(holeStart)~=length(holeEnd)
    holeEnd=[holeEnd,size(flagTemp,2)-1];
end

for ii=1:length(holeStart)
    flagColStart=flagTemp(:,holeStart(ii)-1);
    flagColEnd=flagTemp(:,holeEnd(ii)+1);

    minStart=min(find(flagColStart==7));
    minEnd=min(find(flagColEnd==7));
    
    holeLength=holeEnd(ii)-holeStart(ii);
    newVec=holeStart(ii):holeEnd(ii);
    
    x=[holeStart(ii),holeEnd(ii)];
    v=[minStart,minEnd];
    
    dataFill=round(interp1(x,v,newVec));
    
    for jj=1:length(dataFill)
        flagTemp(dataFill(jj):end,newVec(jj))=7;
    end
end

extInd=find(any(data.FLAG==3,1));
%%%%%%%%%%%%%% Start here
for ii=1:length(extInd)
    dataCol=refl(:,extInd(ii));
    dataMed=median(dataCol,'omitnan');
end
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
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
   % print([figdir,'refl_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

