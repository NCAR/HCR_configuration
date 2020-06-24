% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,4,20,0,0);
% endTime=datetime(2018,2,5,20,0,0);

startTime=datetime(2019,8,7,16,30,0);
endTime=datetime(2019,8,7,17,5,0);

% startTime=datetime(2019,8,11,0,35,0);
% endTime=datetime(2019,8,11,23,55,0);

if etime(datevec(endTime),datevec(startTime))<=1800
    plotPartFlight=1;
else
    plotPartFlight=0;
end

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

figdir=['/h/eol/romatsch/hcrCalib/clouds/brightBand/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data.DBZ=[];
data.LDR=[];
data.VEL_CORR=[];
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

%% Find melting layer
data.dbzMasked=data.DBZ;
data.dbzMasked(data.FLAG>1)=nan;

findMelt=f_meltLayer(data);
zeroInds=find(findMelt==0);
oneInds=find(findMelt==1);
twoInds=find(findMelt==2);
threeInds=find(findMelt==3);

%% Plot

timeMat=repmat(data.time,size(data.TEMP,1),1);

%close all

if plotPartFlight
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);
    
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
    print([figdir,'meltReflShort_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');

%% Plot LDR and VEL

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
    
    ax1=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,data.VEL_CORR,'edgecolor','none');
    view(2);
    caxis([-5 5]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('VEL')
    grid on
        
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'LDR_VEL_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
else
%% Plot longer data stretch
    close all
    
    % Resample for plotting
    newInds=1:100:length(data.time);
    newDBZ=data.DBZ(:,newInds);
    newLDR=data.LDR(:,newInds);
    newVEL=data.VEL_CORR(:,newInds);
    newASL=data.asl(:,newInds);
    newTEMP=data.TEMP(:,newInds);
    newTime=data.time(newInds);
    
    fig3=figure('DefaultAxesFontSize',11,'position',[100,1300,2500,1000]);
    
    s1=subplot(3,1,1);
    hold on;
    sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    colIndsAll=1:length(data.time);
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
        
    s2=subplot(3,1,2);
    hold on;
    sub2=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
    view(2);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    
    s3=subplot(3,1,3);
    hold on;
    sub2=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
    view(2);
    s3.Colormap=jet;
    colorbar
    caxis([-5 5]);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('VEL')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'meltReflLong_',datestr(data.time(1),formatOut),'_to_',datestr(data.time(end),formatOut)],'-dpng','-r0');
end

%% Calculate difference between zero degree and detected melting layer
diffInd=[];

for ii=1:size(findMelt,2)
    if data.elevation(ii)<-85
        bbRay=findMelt(:,ii);
        foundInd=min(find(bbRay==1));
        if ~isempty(foundInd)
            zeroInd=min(find(bbRay==0));
            diffInd=[diffInd foundInd-zeroInd];
        end
    end
end

meanDiff=mean(diffInd);
meanMeters=meanDiff*(data.range(2)-data.range(1));

disp(['Melting layer is on average ',num2str(meanDiff),' gates or ',num2str(meanMeters),' m below the zero degree level.'])